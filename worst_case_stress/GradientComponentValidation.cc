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
#include <LpHoleInflator.hh>
#include <BoundaryPerturbationInflator.hh>

#include <PatternOptimizationJob.hh>
#include <PatternOptimizationConfig.hh>

#include "WCStressOptimizationIterate.hh"
#include "BoundaryPerturbationIterate.hh"

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

void usage(int exitVal, const po::options_description &visibleOptions) {
    cout << "Usage: GradientComponentValidation [options] job.opt component_idx -l'lower_bd' -u'upper_bd' nsamples" << endl;
    cout << visibleOptions << endl;
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

    po::options_description sweepOptions("Sweep Options");
    sweepOptions.add_options()
        ("lower_bd,l",    po::value<double>(),                             "sweep lower bound (must be non-positional to support negative values)")
        ("upper_bd,u",    po::value<double>(),                             "sweep upper bound (must be non-positional to support negative values)")
        ;

    po::options_description patternOptions;
    patternOptions.add_options()
        ("pattern,p",    po::value<string>(),                            "Pattern wire mesh (.obj|wire), or initial mesh for BoundaryPerturbationInflator")
        ("inflator,i",   po::value<string>()->default_value("original"), "Which inflator to use: original (default), lphole, boundary_perturbation")
        ("isotropicParameters,I",                                        "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                            "Use vertex thickness instead of edge thickness (3D only)")
        ("cell_size,c",  po::value<double>()->default_value(5.0),        "Inflation cell size (3D only)")
        ;

    po::options_description meshingOptions("Meshing Options");
    meshingOptions.add_options()
        ("meshingOptions,M", po::value<string>(),                    "Meshing options configuration file")
        ("max_volume,v",     po::value<double>(),                    "Maximum element area for remeshing (overrides meshing options)")
        ("hole_segments",    po::value<size_t>()->default_value(64), "Number of segments in hole boundary (LpHoleInflator)")
        // ("subdivide,S",   po::value<size_t>()->default_value(0),           "Number of subdivisions to run for 3D inflator")
        // ("sub_algorithm,A", po::value<string>()->default_value("simple"), "Subdivision algorithm for 3D inflator (simple or loop)")
        ;

    po::options_description optimizerOptions("Optimizer Options");
    optimizerOptions.add_options()
        ("nIters,n",     po::value<size_t>()->default_value(20),                 "number of iterations")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient step size")
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs")
        // ("vtxNormalPerturbationGradient,N",                                      "use the vertex-normal-based versino of the boundary perturbation gradient")
        ;

    po::options_description objectiveOptions("Objective Options");
    objectiveOptions.add_options()
        ("ignoreShear",                                                   "Ignore the shear components in the isotropic tensor fitting")
        ("pnorm,P",      po::value<double>()->default_value(1.0),         "pnorm used in the Lp global worst case stress measure")
        ("usePthRoot,R",                                                  "Use the true Lp norm for global worst case stress measure (applying pth root)")
        ("WCSWeight",    po::value<double>()->default_value(1.0),         "Weight for the WCS term of the objective")
        ("JSWeight",     po::value<double>()->default_value(0.0),         "Weight for the JS term of the objective")
        ("JVolWeight",   po::value<double>()->default_value(0.0),         "Weight for the JVol term of the objective")
        ("LaplacianRegWeight,r", po::value<double>()->default_value(0.0), "Weight for the boundary Laplacian regularization term")
        ;

    po::options_description elasticityOptions("Elasticity Options");
    elasticityOptions.add_options()
        ("material,m",   po::value<string>(),                    "Base material")
        ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
        ;

    po::options_description generalOptions;
    generalOptions.add_options()
        ("help,h",                                               "Produce this help message")
        ("output,o",     po::value<string>(),                    "Output mesh and fields at each iteration")
        ("volumeMeshOut", po::value<string>(),                   "Output volume mesh at each iteration")
        ("dumpShapeDerivatives"  , po::value<string>(),          "Dump shape derivative fields for JVol, JS, and WCS")
        ;

    po::options_description visibleOptions;
    visibleOptions.add(sweepOptions).add(patternOptions).add(meshingOptions).add(optimizerOptions)
                  .add(objectiveOptions).add(generalOptions).add(elasticityOptions);

    po::options_description cli_opts;
    cli_opts.add(visibleOptions).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visibleOptions);
    }

    bool fail = false;
    if (vm.count("job") + vm.count("component_idx")  + vm.count("lower_bd") + vm.count("upper_bd") + vm.count("nsamples") != 5) {
        cout << "Error: must specify input job.opt file, sweep component index, sweep bounds, and sweep samples" << endl;
        fail = true;
    }

    string inflator = vm["inflator"].as<string>();
    if ((inflator != "lphole") && (vm.count("pattern") == 0)) {
        cout << "Error: must specify pattern mesh" << endl;
        fail = true;
    }
    else if (vm.count("pattern")) {
        std::cerr << "WARNING: pattern argument not expected for LpHoleInflator" << std::endl;
    }

    size_t d = vm["degree"].as<size_t>();
    if (d < 1 || d > 2) {
        cout << "Error: FEM Degree must be 1 or 2" << endl;
        fail = true;
    }

    // set<string> subdivisionAlgorithms = {"simple", "loop"};
    // if (subdivisionAlgorithms.count(vm["sub_algorithm"].as<string>()) == 0) {
    //     cout << "Illegal subdivision algorithm specified" << endl;
    //     fail = true;
    // }

    if (fail || vm.count("help"))
        usage(fail, visibleOptions);

    return vm;
}

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
typedef ScalarField<Real> SField;

#include "InflatorTraits.inl"

template<class _Iterate, class _Inflator, class _Objective>
void genAndReportIterate(_Inflator &inflator, const SField &params, _Objective &fullObjective,
                        size_t compIdx, size_t i, const po::variables_map &args) {
    _Iterate it(inflator, params.domainSize(), &params[0], fullObjective);

    std::cout << i << "\t" << params[compIdx] << "\t"
              << it.evaluateJS() << "\t" << it.gradp_JS()[compIdx] << "\t"
              << it.evaluateWCS() << "\t";
    std::cout << it.gradientWCS_direct_component(compIdx) << "\t";
    std::cout << it.gradientWCS_adjoint()[compIdx]
              << std::endl;

    if (args.count("output"))        it.writeMeshAndFields(    args["output"].as<string>() + "_" + std::to_string(i));
    if (args.count("volumeMeshOut")) it.writeVolumeMesh(args["volumeMeshOut"].as<string>() + "_" + std::to_string(i) + ".msh");

    // TODO: write shape-derivative output as per-element-vertex shape velocity
    // (Should already be supported by MSHFieldWriter)
    // if (args.count("dumpShapeDerivatives")) {
    //     auto path = args["dumpShapeDerivatives"].as<string>()
    //                 + "_" + std::to_string(i) + ".msh";
    //     MSHBoundaryFieldWriter bdryWriter(path, it.mesh());
    //     bdryWriter.addField("dJS",   bpi->boundaryVectorField(it.gradp_JS()),            DomainType::PER_NODE);
    //     bdryWriter.addField("dWCS",  bpi->boundaryVectorField(it.gradientWCS_adjoint()), DomainType::PER_NODE);
    //     bdryWriter.addField("dJVol", bpi->boundaryVectorField(it.gradp_JVol()), DomainType::PER_NODE);
    // }
}

template<size_t _N, size_t _FEMDegree, template<size_t> class _ITraits>
void execute(const po::variables_map &args, const PatternOptimization::Job<_N> *job)
{
    auto inflator_ptr = _ITraits<_N>::construct(args, job);
    auto &inflator = *inflator_ptr;

    // Only relevant to isosurface inflator...
    PatternOptimization::Config::get().useSDNormalShapeVelocityDirectly =
        args.count("directSDNormalVelocity");

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    if (job->numParams() != inflator.numParameters()) {
        for (size_t i = 0; i < inflator.numParameters(); ++i) {
            auto ptype = inflator.parameterType(i);
            cout << "param " << i << " role: " <<
                (ptype == ParameterType::Offset ? "Offset" :
                 ((ptype == ParameterType::Thickness) ? "Thickness" : "Blending"))
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
    size_t compIdx = args["component_idx"].as<size_t>();
    assert(compIdx < params.domainSize());
    double lowerBound = args["lower_bd"].as<double>();
    double upperBound = args["upper_bd"].as<double>();
    size_t nSamples = args["nsamples"].as<size_t>();

    // Create scalarized multi-objective with weights specified by the
    // arguments.
    WCStressOptimization::Objective<_N> fullObjective(targetS,
                                args[  "JSWeight"].as<double>(),
                                args[ "WCSWeight"].as<double>(),
                                args["JVolWeight"].as<double>(),
                                args["LaplacianRegWeight"].as<double>());

    auto &config = WCStressOptimization::Config::get();
    config.globalObjectivePNorm = args["pnorm"].as<double>();

    using LpIterate = WCStressOptimization::Iterate<Simulator>;
    using LinfIterate = WCStressOptimization::Iterate<Simulator,
                            IntegratedWorstCaseObjective<_N, WCStressIntegrandLinf>>;

    for (size_t i = 0; i < nSamples; ++i) {
        params[compIdx] = lowerBound + ((nSamples == 1) ? 0.0
                        : (upperBound - lowerBound) * (double(i) / (nSamples - 1)));
        if (std::isinf(config.globalObjectivePNorm)) genAndReportIterate<LinfIterate>(inflator, params, fullObjective, compIdx, i, args);
        else                                         genAndReportIterate<  LpIterate>(inflator, params, fullObjective, compIdx, i, args);
    }

    BENCHMARK_REPORT();
}

int main(int argc, const char *argv[])
{
    cout << setprecision(20);

    po::variables_map args = parseCmdLine(argc, argv);

    auto job = PatternOptimization::parseJobFile(args["job"].as<string>());

    auto inflator = args["inflator"].as<string>();

    size_t deg = args["degree"].as<size_t>();
    if (auto job2D = dynamic_cast<PatternOptimization::Job<2> *>(job)) {
        if (inflator == "original") {
            if (deg == 1) execute<2, 1, InflatorTraitsConstrainedInflator>(args, job2D);
            if (deg == 2) execute<2, 2, InflatorTraitsConstrainedInflator>(args, job2D);
        }
        else if (inflator == "lphole") {
            if (deg == 1) execute<2, 1, InflatorTraitsLpHole>(args, job2D);
            if (deg == 2) execute<2, 2, InflatorTraitsLpHole>(args, job2D);
        }
        else if (inflator == "boundary_perturbation") {
            if (deg == 1) execute<2, 1, InflatorTraitsBoundaryPerturbation>(args, job2D);
            if (deg == 2) execute<2, 2, InflatorTraitsBoundaryPerturbation>(args, job2D);
        }
        else throw std::runtime_error("Unknown inflator: " + inflator);
    }
    else if (auto job3D = dynamic_cast<PatternOptimization::Job<3> *>(job)) {
        if (inflator == "original") {
            if (deg == 1) execute<3, 1, InflatorTraitsConstrainedInflator>(args, job3D);
            if (deg == 2) execute<3, 2, InflatorTraitsConstrainedInflator>(args, job3D);
        }
        else if (inflator == "lphole") throw std::runtime_error("3D LpHole inflator unimplemented");
        else if (inflator == "boundary_perturbation") {
            if (deg == 1) execute<3, 1, InflatorTraitsBoundaryPerturbation>(args, job3D);
            if (deg == 2) execute<3, 2, InflatorTraitsBoundaryPerturbation>(args, job3D);
        }
        else throw std::runtime_error("Unknown inflator: " + inflator);
    }
    else throw std::runtime_error("Invalid job file.");

    return 0;
}
