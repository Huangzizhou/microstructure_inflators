////////////////////////////////////////////////////////////////////////////////
// PatternOptimization_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Evolves the microstructure parameters to bring the structure's
//      homogenized elasticity tensor closer to a target tensor.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/12/2014 01:15:28
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <PeriodicHomogenization.hh>

#include <WireInflator2D.h>

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "PatternOptimization.hh"
#include "PatternOptimizationJob.hh"

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: PatternOptimization_cli [options] job.opt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job", po::value<string>(), "job configuration file")
        ;
    po::positional_options_description p;
    p.add("job", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help", "Produce this help message")
        ("pattern,p",   po::value<string>(), "Pattern wire mesh (.obj)")
        ("material,m",  po::value<string>(), "base material")
        ("output,o",    po::value<string>()->default_value(""), "output .js mesh + fields at each iteration")
        ("max_area,a",  po::value<double>()->default_value(0.0001), "max_area parameter for wire inflator")
        ("solver",      po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs, levenberg_marquardt")
        ("step,s",      po::value<double>()->default_value(0.0001), "gradient step size")
        ("nIters,n",    po::value<size_t>()->default_value(20), "number of iterations")
        ;

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    bool fail = false;
    if (vm.count("job") == 0) {
        cout << "Error: must specify input job.opt file" << endl;
        fail = true;
    }

    if (vm.count("pattern") == 0) {
        cout << "Error: must specify pattern mesh" << endl;
        fail = true;
    }

    set<string> solvers = {"gradient_descent", "bfgs", "lbfgs", "levenberg_marquardt"};
    if (solvers.count(vm["solver"].as<string>()) == 0) {
        cout << "Illegal solver specified" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

template<size_t _N>
using ETensor = typename LinearElasticityND<_N>::ETensor;
template<size_t _N>
using VField = typename LinearElasticityND<_N>::VField;
typedef ScalarField<Real> SField;

template<size_t _N>
void execute(const po::variables_map &args,
             const Job<_N> *job);

template<>
void execute<2>(const po::variables_map &args,
                const Job<2> *job)
{
	WireInflator2D inflator(args["pattern"].as<string>());
    TessellationParameters t_params;
    t_params.max_area = args["max_area"].as<double>();

    constexpr size_t _N = 2;
    auto targetC = job->targetMaterial.getTensor();
    ETensor<2> targetS = targetC.inverse();

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    size_t nParams = inflator.patternGenerator().numberOfParameters();
    if (job->numParams() != nParams)
        throw runtime_error("Invalid number of parameters.");

    // Set up simulator's (base) material
    auto &mat = LinearElasticityND<_N>::
        template homogenousMaterial<Materials::Constant>();
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    SField params(nParams);
    params = job->initialParams;

    Optimizer<_N> optimizer(inflator, t_params, job->radiusBounds,
            job->translationBounds);
    string solver = args["solver"].as<string>(),
           output = args["output"].as<string>();
    size_t niters = args["nIters"].as<size_t>();
    if (solver == "levenberg_marquardt")
        optimizer.optimize_lm(params, targetS, output);
    else if (solver == "gradient_descent")
        optimizer.optimize_gd(params, targetS, niters,
                          args["step"].as<double>(), output);
    else if (solver == "bfgs")
        optimizer.optimize_bfgs(params, targetS, niters, output);
    else if (solver == "lbfgs")
        optimizer.optimize_bfgs(params, targetS, niters, output, 10);
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    cout << setprecision(16);

    po::variables_map args = parseCmdLine(argc, argv);

    auto job = parseJobFile(args["job"].as<string>());

    if (auto job2D = dynamic_cast<Job<2> *>(job)) {
        execute(args, job2D);
    }
    // else if (auto job3D = dynamic_cast<Job<3> *>(job)) {
    //  execute(args, job3D);
    // }
    else throw std::runtime_error("Only 2D jobs currently supported.");

    return 0;
}
