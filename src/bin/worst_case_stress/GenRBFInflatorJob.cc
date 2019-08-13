////////////////////////////////////////////////////////////////////////////////
// GenRBFInflatorJob.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Create a job for the RBF inflator.
*/
////////////////////////////////////////////////////////////////////////////////

#include <MeshFEM/Future.hh>
#include <MeshFEM/StringUtils.hh>
#include <inflators/wrappers/RBFInflator.hh>
#include <inflators/wrappers/RBFOrthoInflator.hh>
#include <pattern_optimization/PatternOptimizationJob.hh>

#include <iostream>
#include <cstdlib>
#include <string>
#include <cassert>
#include <stdexcept>
#include <cmath>

#include <CLI/CLI.hpp>

using namespace std;

struct Args {
    size_t dim;
    double maxCoeff;
    double targetVolume;
    std::string elasticityTensor;
    std::string targetPng;
    std::string initialPng;
    bool ortho = false;
};

int parseCmdLine(int argc, const char *argv[], Args &args)
{
    // Parse arguments
    CLI::App app{"GenRBFInflatorJob"};

    app.add_option("dim",                args.dim,               "dimension");
    app.add_option("--maxCoeff",         args.maxCoeff,          "maximum absolute value of a coefficient");
    app.add_option("--targetVolume",     args.targetVolume,      "target volume");
    app.add_option("--elasticityTensor", args.elasticityTensor,  "target tensor specifier (Young,Poisson)");
    app.add_option("--targetPng",        args.targetPng,         "target shape");
    app.add_option("--initialPng",       args.initialPng,        "initial shape");
    app.add_flag(  "--ortho",            args.ortho,             "orthotropic symmetry (uses less parameters). Set initial/target parameters accordingly");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    return 0;
}

vector<Real> parseVecArg(const string &arg, size_t N) {
    vector<string> components;

    components = MeshFEM::split(arg, ",");
    //boost::split(components, arg, boost::is_any_of(","));
    if (components.size() != N) throw runtime_error("Invalid number of components in argument " + arg);

    vector<Real> result;
    for (const string &c : components)
        result.push_back(stod(c));
    return result;
}

int main(int argc, const char *argv[])
{
    Args args;
    args.elasticityTensor = "1,0";
    args.maxCoeff = 10.0;
    int exitCode = parseCmdLine(argc, argv, args);
    if (exitCode) {
        return exitCode;
    }

    vector<Real> targetModuli = parseVecArg(args.elasticityTensor, 2);

    unique_ptr<PatternOptimization::JobBase> job;
    auto job2D = Future::make_unique<PatternOptimization::Job<2>>();
    job2D->targetMaterial.setIsotropic(targetModuli[0], targetModuli[1]);
    job = move(job2D);

    job->numberCustomTypes = 1;

    size_t dim = args.dim;

    if (args.maxCoeff > 0.0) {
        double maxCoeff = args.maxCoeff;

        for (size_t p = 0; p < dim*dim; ++p) {
            job->varLowerBounds.emplace(p, -maxCoeff);
            job->varUpperBounds.emplace(p,  maxCoeff);
        }

        job->custom1Bounds = {-maxCoeff, maxCoeff};
    }

    // Create inflator
    if (!args.targetPng.empty()) {
        std::unique_ptr<Inflator<2>> inflator;
        if (args.ortho) {
            size_t d = args.dim;
            Real epsilon = (d + d - 1.0) / 2.0;
            inflator = Future::make_unique<RBFOrthoInflator>(args.targetPng, epsilon, d);
        } else {
            size_t d = args.dim;
            Real epsilon = d / 2;
            inflator = Future::make_unique<RBFInflator>(args.targetPng, epsilon, d);
        }

        job->targetParams = inflator->defaultParameters();
    }

    if (!args.initialPng.empty()) {
        std::unique_ptr<Inflator<2>> inflator;
        if (args.ortho) {
            size_t d = args.dim;
            Real epsilon = (d + d - 1.0) / 2.0;
            inflator = Future::make_unique<RBFOrthoInflator>(args.initialPng, epsilon, d);
        } else {
            size_t d = args.dim;
            Real epsilon = d / 2;
            inflator = Future::make_unique<RBFInflator>(args.initialPng, epsilon, d);
        }

        job->initialParams = inflator->defaultParameters();
    }


    if (args.targetVolume) {
        job->targetVolume = args.targetVolume;
    }

    // Fake value for other bounds
    job->translationBounds = {0.1, 0.8};
    job->radiusBounds   = {0.1, 0.2};
    job->blendingBounds = {0.0, 0.1};

    job->writeJobFile(cout);

    return 0;
}
