////////////////////////////////////////////////////////////////////////////////
// RandomJob.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Generates a random job trying to reach a random random point in the feasible
//  region from another random point.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/16/2014 14:23:18
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
    cout << "Usage: RandomJob [options] config.opt out_job.opt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job",     po::value<string>(), "job configuration file")
        ("out_job", po::value<string>(), "out job configuration file")
        ;
    po::positional_options_description p;
    p.add("job", 1);
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

    if (vm.count("out_job") == 0) {
        cout << "Error: must specify output file out_job.opt" << endl;
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

    double md = vm["max_distance"].as<double>();
    if (md < 0.0) {
        cout << "Error: maximum distance must be creater than 0." << endl;
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

    // Generate random ground-truth pattern parameters
    Job<_N> outJob(*job);
    outJob.trueParams.resize(inflator.numParameters());

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> radiusParameter(job->radiusBounds[0], job->radiusBounds[1]),
                                     translationParameter(job->translationBounds[0], job->translationBounds[1]);
    for (size_t i = 0; i < inflator.numParameters(); ++i) {
        if (inflator.parameterType(i) == ParameterType::Thickness)
            outJob.trueParams[i] = radiusParameter(gen);
        else 
            outJob.trueParams[i] = translationParameter(gen);
    }
    
    // Find corresponding elasticity tensor
    Iterate<Simulator> randomIter(inflator, outJob.trueParams.size(),
                                  &outJob.trueParams[0], ETensor<_N>());
    outJob.targetMaterial.setTensor(randomIter.elasticityTensor());

    // Generate random initial pattern parameters within a ball of specified
    // radius around the ground truth point in pattern parameter space.
    // (But only pick feasible points).
    // To speed up search, we hone in on a tight box containing all allowed points.
    Real translationRange = job->translationBounds[1] - job->translationBounds[0];
    Real radiusRange      =      job->radiusBounds[1] -      job->radiusBounds[0];
    Real maxDistance = args["max_distance"].as<double>();

    // Keep picking from this box until we get a valid point.
    uniform_real_distribution<> relativeOffsetComponent(std::max(-maxDistance, -1.0),
                                                        std::min(maxDistance,   1.0));
    Eigen::VectorXd offset(inflator.numParameters());
    outJob.initialParams.resize(inflator.numParameters());
    do {
        for (size_t i = 0; i < inflator.numParameters(); ++i) {
            Real param;
            if (inflator.parameterType(i) == ParameterType::Thickness) {
                do {
                    offset(i) = relativeOffsetComponent(gen);
                    param = outJob.trueParams[i] + radiusRange * offset(i);
                } while ((param < job->radiusBounds[0]) || (param > job->radiusBounds[1]));
            }
            else {
                do {
                    offset(i) = relativeOffsetComponent(gen);
                    param = outJob.trueParams[i] + translationRange * offset(i);
                } while ((param < job->translationBounds[0]) || (param > job->translationBounds[1]));
            }
            outJob.initialParams[i] = param;
        }
    } while (offset.norm() > maxDistance);

    outJob.writeJobFile(args["out_job"].as<string>());

    BENCHMARK_REPORT();
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    cout << setprecision(30);

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
