////////////////////////////////////////////////////////////////////////////////
// GradientComponentValidation.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Validates a single component of the pattern optimization gradient. Does
//      this by doing a fine 1D sweep of the objective around the sample point,
//      writing the objective and gradient component for each sample.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/26/2014 14:57:55
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include <GlobalBenchmark.hh>

#include <Inflator.hh>
#include <PatternOptimizationJob.hh>
#include <PatternOptimizationConfig.hh>
#include "WCStressOptimizationIterate.hh"

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

namespace po = boost::program_options;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: GradientComponentValidation [options] job.opt component_idx -l'lower_bd' -u'upper_bd' nsamples" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job", po::value<string>(), "job configuration file")
        ("component_idx", po::value<size_t>(), "index of component to sweep")
        ("nsamples", po::value<size_t>(), "number of steps to sweep")
        ;
    po::positional_options_description p;
    p.add("job", 1);
    p.add("component_idx", 1);
    p.add("nsamples", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help",        "Produce this help message")
        ("lower_bd,l", po::value<double>(), "sweep lower bound (must be positional to support negative values)")
        ("upper_bd,u", po::value<double>(), "sweep upper bound (must be positional to support negative values)")
        ("pattern,p",       po::value<string>(), "Pattern wire mesh (.obj|wire)")
        ("material,m",      po::value<string>(), "base material")
        ("degree,d",        po::value<size_t>()->default_value(2),        "FEM Degree")
        ("output,o",        po::value<string>()->default_value(""),       "output .js mesh + fields at each iteration")
        ("matrixOut,O",     po::value<string>()->default_value(""),       "output matrix at each iteration")
        ("volumeMeshOut",   po::value<string>()->default_value(""),       "output volume mesh at each iteration")
        ("dofOut",          po::value<string>()->default_value(""),       "output pattern dofs in James' format at each iteration")
        ("subdivide,S",  po::value<size_t>()->default_value(0),           "number of subdivisions to run for 3D inflator")
        ("sub_algorithm,A", po::value<string>()->default_value("simple"), "subdivision algorithm for 3D inflator (simple or loop)")
        ("max_volume,v", po::value<double>(),                             "maximum element volume parameter for wire inflator")
        ("cell_size,c",  po::value<double>()->default_value(5.0),         "Inflation cell size (3D only)")
        ("directSDNormalVelocity,n",                                      "For isosurface inflator: use signed distance normal velocity directly (instead of multiplying by vertex normal . face normal)")
        ("isotropicParameters,I",                                         "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                             "Use vertex thickness instead of edge thickness (3D only)")
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
    if (vm.count("job") + vm.count("component_idx")  + vm.count("lower_bd") + vm.count("upper_bd") + vm.count("nsamples") != 5) {
        cout << "Error: must specify input job.opt file, sweep component index, sweep bounds, and sweep samples" << endl;
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
typedef ScalarField<Real> SField;

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args, const PatternOptimization::Job<_N> *job)
{
    ConstrainedInflator<_N> inflator(job->parameterConstraints,
                args["pattern"].as<string>(),
                args["cell_size"].as<double>(),
                0.5 * sqrt(2),
                args.count("isotropicParameters"),
                args.count("vertexThickness"));
    if (args.count("max_volume"))
        inflator.setMaxElementVolume(args["max_volume"].as<double>());
    if (_N == 3) {
        inflator.configureSubdivision(args["sub_algorithm"].as<string>(),
                                      args["subdivide"].as<size_t>());
    }

    // Only relevant to isosurface inflator...
    PatternOptimization::Config::get().useSDNormalShapeVelocityDirectly =
        args.count("directSDNormalVelocity");

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    // cout << "target tensor: " << targetC << endl;

    if (job->numParams() != inflator.numParameters()) {
        for (size_t i = 0; i < inflator.numParameters(); ++i) {
            cout << "param " << i << " role: " <<
                (inflator.parameterType(i) == ParameterType::Offset ? "Offset" : "Thickness")
                << endl;
        }
        throw runtime_error("Invalid number of parameters.");
    }

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    SField params(job->initialParams);
    string output = args["output"].as<string>();
    string matrixOutput = args["matrixOut"].as<string>();
    string volumeMeshOut = args["volumeMeshOut"].as<string>();
    size_t compIdx = args["component_idx"].as<size_t>();
    assert(compIdx < params.domainSize());
    double lowerBound = args["lower_bd"].as<double>();
    double upperBound = args["upper_bd"].as<double>();
    size_t nSamples = args["nsamples"].as<size_t>();

    string dofOut = args["dofOut"].as<string>();
    if (dofOut != "")
        inflator.setDoFOutputPrefix(dofOut);

    for (size_t i = 0; i < nSamples; ++i) {
        params[compIdx] = lowerBound + ((nSamples == 1) ? 0.0
                        : (upperBound - lowerBound) * (double(i) / (nSamples - 1)));
        WCStressOptimization::Iterate<Simulator> it(inflator, params.domainSize(), &params[0], targetS);
        std::cout << i << "\t" << params[compIdx] << "\t"
                  << it.evaluateJS() << "\t" << it.gradp_JS()[compIdx] << "\t"
                  << it.evaluateWCS() << "\t" << it.gradientWCS_direct()[compIdx] << "\t"
                  << it.gradientWCS_adjoint()[compIdx]
                  << std::endl;
        if (output != "")          it.writeMeshAndFields(       output + "_" + std::to_string(i));
        if (matrixOutput != "")  it.dumpSimulationMatrix( matrixOutput + "_" + std::to_string(i) + ".txt");
        if (volumeMeshOut != "")      it.writeVolumeMesh(volumeMeshOut + "_" + std::to_string(i) + ".msh");
    }

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
    cout << setprecision(20);

    po::variables_map args = parseCmdLine(argc, argv);

    auto job = PatternOptimization::parseJobFile(args["job"].as<string>());

    size_t deg = args["degree"].as<size_t>();
    if (auto job2D = dynamic_cast<PatternOptimization::Job<2> *>(job)) {
        if (deg == 1) execute<2, 1>(args, job2D);
        if (deg == 2) execute<2, 2>(args, job2D);
    }
    else if (auto job3D = dynamic_cast<PatternOptimization::Job<3> *>(job)) {
        if (deg == 1) execute<3, 1>(args, job3D);
        if (deg == 2) execute<3, 2>(args, job3D);
    }
    else throw std::runtime_error("Invalid job file.");

    return 0;
}
