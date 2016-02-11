////////////////////////////////////////////////////////////////////////////////
// WCSOptimization_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Evolves the microstructure parameters to bring the structure's
//      homogenized elasticity tensor closer to a target while also minimizing
//      worst-case stress.
//
//      Example invocation using BoundaryPerturbationInflator:
//      ./WCSOptimization_cli -p demo.msh
//              -m $MICRO_DIR/materials/B9Creator.material ./test_2D_job.opt
//              --alpha 0 --inflator boundary_perturbation -o iterate -d1
//              --step 1e-21 -P4 -n 80
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/12/2014 01:15:28
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include <GlobalBenchmark.hh>

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

#include <Inflator.hh>
#include <PatternOptimizationJob.hh>
#include <PatternOptimizationConfig.hh>

#include "WCSOptimization.hh"
#include "WCSObjective.hh"

#include "WCStressOptimizationIterate.hh"
#include "BoundaryPerturbationIterate.hh"

#include <LpHoleInflator.hh>
#include <BoundaryPerturbationInflator.hh>

namespace po = boost::program_options;
using namespace std;
using namespace WCStressOptimization;

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

    po::options_description patternOptions;
    patternOptions.add_options()
        ("pattern,p",    po::value<string>(),                            "Pattern wire mesh (.obj|wire), or initial mesh for BoundaryPerturbationInflator")
        ("inflator,i",   po::value<string>()->default_value("original"), "Which inflator to use: original (default), lphole, boundary_perturbation")
        ("isotropicParameters,I",                                        "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                            "Use vertex thickness instead of edge thickness (3D only)")
        ("cell_size,c",  po::value<double>()->default_value(5.0),        "Inflation cell size (3D only)")
        ;

    po::options_description meshingOptions;
    meshingOptions.add_options()
        ("max_volume,v",   po::value<double>(), "Maximum element area for remeshing (overrides meshing options)")
        ;

    po::options_description optimizerOptions;
    optimizerOptions.add_options()
        ("nIters,n",     po::value<size_t>()->default_value(20),                 "number of iterations")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient step size")
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs")
        ("vtxNormalPerturbationGradient,N",                                      "use the vertex-normal-based versino of the boundary perturbation gradient")
        ;

    po::options_description objectiveOptions;
    objectiveOptions.add_options()
        ("ignoreShear",                                                   "Ignore the shear components in the isotropic tensor fitting")
        ("pnorm,P",      po::value<double>()->default_value(1.0),         "pnorm used in the Lp global worst case stress measure")
        ("usePthRoot,R",                                                  "use the true Lp norm for global worst case stress measure (applying pth root)")
        ("WCSWeight",    po::value<double>()->default_value(1.0),         "Weight for the WCS term of the objective")
        ("JSWeight",     po::value<double>()->default_value(0.0),         "Weight for the JS term of the objective")
        ("JVolWeight",   po::value<double>()->default_value(0.0),         "Weight for the JVol term of the objective")
        ("LaplacianRegWeight,r", po::value<double>()->default_value(0.0), "Weight for the boundary Laplacian regularization term")
        ;

    po::options_description generalOptions;
    generalOptions.add_options()
        ("help,h",                                               "Produce this help message")
        ("material,m",   po::value<string>(),                    "base material")
        ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
        ("output,o",     po::value<string>()->default_value(""), "output .js mesh + fields at each iteration")
        ("dumpShapeDerivatives"  , po::value<string>(),          "Dump shape derivative fields for JVol, JS, and WCS")
        ;

    po::options_description visibleOptions;
    visibleOptions.add(patternOptions).add(meshingOptions).add(optimizerOptions)
                  .add(objectiveOptions).add(generalOptions);

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
    if (vm.count("job") == 0) {
        cout << "Error: must specify input job.opt file" << endl;
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

    set<string> solvers = {"gradient_descent", "bfgs", "lbfgs"};
    if (solvers.count(vm["solver"].as<string>()) == 0) {
        cout << "Illegal solver specified" << endl;
        fail = true;
    }

    set<string> inflators = {"original", "lphole", "boundary_perturbation"};
    if (inflators.count(vm["inflator"].as<string>()) == 0) {
        cout << "Illegal inflator specified" << endl;
        fail = true;
    }

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

template<size_t _N, size_t _FEMDegree, template<size_t> class _ITraits>
void execute(const po::variables_map &args, const PatternOptimization::Job<_N> *job)
{
    using _Inflator = typename _ITraits<_N>::type;
    auto inflator_ptr = _ITraits<_N>::construct(args, job);

    auto &inflator = *inflator_ptr;

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    cout << "target tensor: " << targetC << endl;

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    // Configure WCS Objective
    auto &wcsConfig = WCStressOptimization::Config::get();
    wcsConfig.globalObjectivePNorm = args["pnorm"].as<double>();
    if (args.count("usePthRoot"))
        wcsConfig.globalObjectiveRoot = 2.0 * wcsConfig.globalObjectivePNorm;
    wcsConfig.useVtxNormalPerturbationGradientVersion = args.count("vtxNormalPerturbationGradient");

    SField params = _ITraits<_N>::initParams(inflator_ptr, args, job);

    PatternOptimization::Config::get().ignoreShear = args.count("ignoreShear");
    if (PatternOptimization::Config::get().ignoreShear) cout << "Ignoring shear components" << endl;
    Optimizer<Simulator, _Inflator, _ITraits<_N>::template Iterate>
        optimizer(inflator, job->radiusBounds,   job->translationBounds, job->blendingBounds,
                            job->varLowerBounds, job->varUpperBounds);

    // Create scalarized multi-objective with weights specified by the
    // arguments.
    WCStressOptimization::Objective<_N> fullObjective(targetS,
                                args[  "JSWeight"].as<double>(),
                                args[ "WCSWeight"].as<double>(),
                                args["JVolWeight"].as<double>(),
                                args["LaplacianRegWeight"].as<double>());

    string solver = args["solver"].as<string>(),
           output = args["output"].as<string>();
    size_t niters = args["nIters"].as<size_t>();

    if (solver == "gradient_descent") optimizer.optimize_gd(  params, fullObjective, niters, args["step"].as<double>(), output);
    else if (solver == "bfgs")        optimizer.optimize_bfgs(params, fullObjective, niters, output, 0);
    else if (solver == "lbfgs")       optimizer.optimize_bfgs(params, fullObjective, niters, output, 10);

    std::vector<Real> result(params.domainSize());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = params[i];

    _ITraits<_N>::finalize(inflator_ptr, result, args, job);

    BENCHMARK_REPORT();
}

int main(int argc, const char *argv[])
{
    cout << setprecision(16);

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
