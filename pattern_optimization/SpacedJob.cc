////////////////////////////////////////////////////////////////////////////////
// SpacedJob.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Generates jobs spaced evenly along a particular parameter. All other
//      parameters are set to the middle of their bound range.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/28/2014 17:07:35
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include "LinearElasticity.hh"
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include "GlobalBenchmark.hh"

#include "Inflator.hh"

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <random>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "PatternOptimizationJob.hh"
#include "PatternOptimizationIterate.hh"

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: SpacedJob [options] config.opt component_idx njobs out_job_prefix" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job",     po::value<string>(), "job configuration file")
        ("component_idx", po::value<size_t>(), "index of component to sweep")
        ("njobs", po::value<size_t>(), "number of jobs to space")
        ("out_job", po::value<string>(), "out job configuration file (jobs out_job_#.opt is written)")
        ;
    po::positional_options_description p;
    p.add("job", 1);
    p.add("component_idx", 1);
    p.add("njobs", 1);
    p.add("out_job", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help",        "Produce this help message")
        ("pattern,p",    po::value<string>(), "Pattern wire mesh (.obj|wire)")
        ("material,m",   po::value<string>(), "base material")
        ("degree,d",     po::value<size_t>()->default_value(2),           "FEM Degree")
        ("subdivide,S",  po::value<size_t>()->default_value(0),           "number of subdivisions to run for 3D inflator")
        ("sub_algorithm,A", po::value<string>()->default_value("simple"), "subdivision algorithm for 3D inflator (simple or loop)")
        ("max_volume,v", po::value<double>(),                             "maximum element volume parameter for wire inflator")
        ("max_distance,D", po::value<double>()->default_value(1),         "maximum distance to true pattern parameters from random initial point (in relative units)")
        ("initAtUpper,U",                                                 "switch to using upper bound (instead of lower) as the initial parameter point")
        ("lowerBd,-l", po::value<double>(), "lower bound of range to subdivide (defaults to lower bd from config.opt)")
        ("upperBd,-u", po::value<double>(), "upper bound of range to subdivide (defaults to upper bd from config.opt)")
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
    if (vm.count("out_job") == 0) {
        cout << "Error: must specify all positional arguments" << endl;
        fail = true;
    }

    if (vm.count("pattern") == 0) {
        cout << "Error: must specify pattern mesh" << endl;
        fail = true;
    }

    size_t d = vm["degree"].as<size_t>();
    if (d < 1 || d > 2) {
        cout << "Error: FEM Degree must be 1 or 2" << endl;
        fail = true;
    }

    set<string> subdivisionAlgorithms = {"simple", "loop"};
    if (subdivisionAlgorithms.count(vm["sub_algorithm"].as<string>()) == 0) {
        cout << "Illegal subdivision algorithm specified" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
template<size_t _N>
using VField = VectorField<Real, _N>;
typedef ScalarField<Real> SField;

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args, const Job<_N> *job)
{
	Inflator<_N> inflator(args["pattern"].as<string>());
    if (args.count("max_volume"))
        inflator.setMaxElementVolume(args["max_volume"].as<double>());
    if (_N == 3) {
        inflator.configureSubdivision(args["sub_algorithm"].as<string>(),
                                      args["subdivide"].as<size_t>());
    }

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());
    
    size_t nParams = inflator.numParameters();
    size_t compIdx = args["component_idx"].as<size_t>();
    if (compIdx >= nParams) throw std::runtime_error("Invalid component_idx");

    // Set ground truth/initial parameters all to midpoint initially
    Job<_N> outJob(*job);
    outJob.trueParams.resize(nParams);
    outJob.initialParams.resize(nParams);
    Real translationMid = (job->translationBounds[0] + job->translationBounds[1]) / 2;
    Real      radiusMid = (     job->radiusBounds[0] +      job->radiusBounds[1]) / 2;

    for (size_t p = 0; p < inflator.numParameters(); ++p) {
        if (inflator.parameterType(p) == ParameterType::Thickness)
            outJob.initialParams[p] = outJob.trueParams[p] = radiusMid;
        else 
            outJob.initialParams[p] = outJob.trueParams[p] = translationMid;
    }

    // Default to using the parameter's constraints as our sweep upper and
    // lower bound, but the {lower,upper}Bd arguments override them.
    double lowerBound, upperBound;
    if (inflator.parameterType(compIdx) == ParameterType::Thickness) { lowerBound =      job->radiusBounds[0]; upperBound =      job->radiusBounds[1]; }
    else                                                             { lowerBound = job->translationBounds[0]; upperBound = job->translationBounds[1]; }
    if (args.count("lowerBd")) lowerBound = args["lowerBd"].as<double>();
    if (args.count("upperBd")) upperBound = args["upperBd"].as<double>();

    size_t nJobs = args["njobs"].as<size_t>();
    double range = upperBound - lowerBound;
    double interval = range / nJobs;
    for (size_t j = 0; j < nJobs; ++j) {
        // True point is at the middle of the interval
        outJob.trueParams[compIdx] = lowerBound + interval * (j + 0.5);
        if (args.count("initAtUpper")) outJob.initialParams[compIdx] = lowerBound + interval * (j + 1.0);
        else                           outJob.initialParams[compIdx] = lowerBound + interval * (j);
        // Find corresponding elasticity tensor
        Iterate<Simulator> iter(inflator, outJob.trueParams.size(),
                                      &outJob.trueParams[0], ETensor<_N>());
        outJob.targetMaterial.setTensor(iter.elasticityTensor());
        outJob.writeJobFile(args["out_job"].as<string>() + std::to_string(j) + ".opt");
    }
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

    size_t deg = args["degree"].as<size_t>();
    if (auto job2D = dynamic_cast<Job<2> *>(job)) {
        if (deg == 1) execute<2, 1>(args, job2D);
        if (deg == 2) execute<2, 2>(args, job2D);
    }
    else if (auto job3D = dynamic_cast<Job<3> *>(job)) {
        if (deg == 1) execute<3, 1>(args, job3D);
        if (deg == 2) execute<3, 2>(args, job3D);
    }
    else throw std::runtime_error("Invalid job file.");

    return 0;
}
