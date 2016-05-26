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
#include <functional>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include <Inflator.hh>
#include <MakeInflator.hh>
#include <PatternOptimizationJob.hh>
#include <PatternOptimizationConfig.hh>

#include <BoundConstraints.hh>
#include <IterateFactory.hh>
#include <IterateManager.hh>

#include <optimizers/ceres.hh>
#include <optimizers/dlib.hh>
#include <optimizers/gradient_descent.hh>

#include <PatternOptimizationIterate.hh>

#include <OptimizerConfig.hh>
#include <PatternOptimizationConfig.hh>
#include <objective_terms/TensorFit.hh>
#include <objective_terms/ProximityRegularization.hh>
#include "WCSObjectiveTerm.hh"

namespace po = boost::program_options;
namespace PO = PatternOptimization;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: WCSOptimization_cli [options] job.opt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

using OptimizerMap =
    map<string, std::function<void(ScalarField<Real> &, const PO::BoundConstraints &,
            PO::IterateManagerBase &, const PO::OptimizerConfig &, const string &)>>;

OptimizerMap optimizers = {
    {"levenberg_marquardt",  optimize_ceres_lm},
    {"dogleg",               optimize_ceres_dogleg},
    {"bfgs",                 optimize_dlib_bfgs},
    {"lbfgs",                optimize_dlib_bfgs},
    {"gradient_descent",     optimize_gd}
};

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
        ("inflator,i",   po::value<string>()->default_value("isosurface"), "Which inflator to use: isosurface (default), lphole, boundary_perturbation")
        ("isotropicParameters,I",                                        "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                            "Use vertex thickness instead of edge thickness (3D only)")
        ("cell_size,c",  po::value<double>(),                            "Inflation cell size (3D only)")
        ;

    po::options_description meshingOptions;
    meshingOptions.add_options()
        ("meshingOptions,M", po::value<string>(),  "Meshing options configuration file")
        ("max_volume,v",     po::value<double>(),  "Maximum element area for remeshing (overrides meshing options)")
        ("hole_segments",    po::value<size_t>(),  "Number of segments in hole boundary for LpHoleInflator (default: 64)")
        ("subdivide,S",  po::value<size_t>(),      "Number of subdivisions to run for James' inflator (default: 0)")
        ("sub_algorithm,A", po::value<string>(),   "Subdivision algorithm for James' inflator (simple or loop, default: simple)")
        ("fullCellInflator",                       "Use the full period cell inflator instead of the reflection-based one")
        ;

    po::options_description optimizerOptions;
    optimizerOptions.add_options()
        ("nIters,n",     po::value<size_t>()->default_value(20),                 "number of iterations")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient step size")
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs")
        ;

    po::options_description objectiveOptions;
    objectiveOptions.add_options()
        ("ignoreShear",                                                   "Ignore the shear components in the isotropic tensor fitting")
        ("pnorm,P",      po::value<double>()->default_value(1.0),         "pnorm used in the Lp global worst case stress measure")
        ("usePthRoot,R",                                                  "Use the true Lp norm for global worst case stress measure (applying pth root)")
        ("WCSWeight",    po::value<double>()->default_value(1.0),         "Weight for the WCS term of the objective")
        ("JSWeight",     po::value<double>(),                             "Use the NLLS tensor fitting term with specified weight.")
        ("proximityRegularizationWeight", po::value<double>(),            "Use a quadratic proximity regularization term with the specified weight.")
        ("LaplacianRegWeight,r", po::value<double>()->default_value(0.0), "Weight for the boundary Laplacian regularization term")
        ;

    po::options_description elasticityOptions;
    elasticityOptions.add_options()
        ("material,m",   po::value<string>(),                    "Base material")
        ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
        ;

    po::options_description gvOptions;
    gvOptions.add_options()
        ("validateGradientComponent", po::value<size_t>(),                   "Run gradient component validation instead of optimization")
        ("nsamples",                  po::value<size_t>()->default_value(5), "Number of gradient component validation samples")
        ("range",                     po::value<string>(),                   "Absolute sweep range (lower:upper)")
        ("rangeRelative",             po::value<double>(),                   "Relative sweep range: current +/- arg * paramBound(compIdx)")
        ("singleIteration",           po::value<size_t>(),                   "Only run the ith iteration of the validation")
        ;

    po::options_description generalOptions;
    generalOptions.add_options()
        ("help,h",                                      "Produce this help message")
        ("output,o",     po::value<string>(),           "output mesh and fields at each iteration")
        ("dumpShapeDerivatives"  , po::value<string>(), "Dump shape derivative fields for JVol, JS, and WCS")
        ;

    po::options_description visibleOptions;
    visibleOptions.add(patternOptions).add(meshingOptions).add(optimizerOptions)
                  .add(objectiveOptions).add(elasticityOptions).add(gvOptions).add(generalOptions);

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

    if (fail || vm.count("help"))
        usage(fail, visibleOptions);

    return vm;
}

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
typedef ScalarField<Real> SField;

template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args, const PO::Job<_N> *job)
{
    auto infl_ptr = make_inflator<_N>(args["inflator"].as<string>(),
                                     filterInflatorOptions(args),
                                     job->parameterConstraints);
    auto &inflator = *infl_ptr;

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();
    bool gradientValidationMode = args.count("validateGradientComponent");
    if (!gradientValidationMode) {
        cout << "Target moduli:\t";
        targetC.printOrthotropic(cout);
        cout << endl;
    }

    // cout << "target tensor: " << targetC << endl;

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    SField params = job->validatedInitialParams(inflator);

    // TODO: Laplacian regularization term (probably only needed for boundary
    // perturbation version.

    using WCSTermConfig       = PO::ObjectiveTerms::IFConfigWorstCaseStress<Simulator>;
    using TensorFitTermConfig = PO::ObjectiveTerms::IFConfigTensorFit<Simulator>;
    using PRegTermConfig      = PO::ObjectiveTerms::IFConfigProximityRegularization;

    // TODO: what to do for boundary perturbation inflator? higher order optimizers should
    // probably check that they're running on a parametric iterate factory, and
    // gradient descent should allow resizing of params otherwise.
    auto ifactory = PO::make_iterate_factory<PO::Iterate<Simulator>,
                                             WCSTermConfig,
                                             TensorFitTermConfig,
                                             PRegTermConfig>(inflator);

    ////////////////////////////////////////////////////////////////////////////
    // Configure the objective terms
    ////////////////////////////////////////////////////////////////////////////
    ifactory->WCSTermConfig      ::enabled = true;
    ifactory->TensorFitTermConfig::enabled = args.count("JSWeight");
    ifactory->PRegTermConfig     ::enabled = args.count("proximityRegularizationWeight");

    // Configure WCS Objective
    // By default, an "Lp norm" objective is really the p^th power of the Lp norm.
    // To use the true "Lp norm", globalObjectiveRoot must be set to
    // 2.0 * globalObjectivePNorm (since pointwise WCS is already squared (e.g. Frobenius) norm)
    ifactory->WCSTermConfig::weight = args["WCSWeight"].as<double>();
    Real pnorm = args["pnorm"].as<double>();
    ifactory->WCSTermConfig::globalObjectivePNorm = pnorm;
    ifactory->WCSTermConfig::globalObjectiveRoot  = args.count("usePthRoot") ? 2.0 * pnorm : 1.0;

    if (args.count("JSWeight")) {
        ifactory->TensorFitTermConfig::weight  = args["JSWeight"].as<double>();
        ifactory->TensorFitTermConfig::targetS = targetS;
        ifactory->TensorFitTermConfig::ignoreShear = args.count("ignoreShear");
        if (ifactory->TensorFitTermConfig::ignoreShear) cout << "Ignoring shear components" << endl;
    }

    if (args.count("proximityRegularizationWeight")) {
        ifactory->PRegTermConfig::initParams = job->initialParams;
        ifactory->PRegTermConfig::weight     = args["proximityRegularizationWeight"].as<double>();
    }

    auto imanager = make_iterate_manager(std::move(ifactory));
    PO::BoundConstraints bdcs(inflator, job->radiusBounds, job->translationBounds, job->blendingBounds,
                              job->varLowerBounds, job->varUpperBounds);

    ////////////////////////////////////////////////////////////////////////////
    // Gradient component validation, if requested, bypasses optimization 
    ////////////////////////////////////////////////////////////////////////////
    if (gradientValidationMode) {
        size_t compIdx = args["validateGradientComponent"].as<size_t>();
        if (compIdx >= params.domainSize()) throw runtime_error("Gradient component index out of bounds");
        if (args.count("range") == args.count("rangeRelative"))
            throw runtime_error("Either range or rangeRelative must be specified (not both)");

        if (!bdcs.hasLowerBound.at(compIdx) || !bdcs.hasUpperBound.at(compIdx))
            throw runtime_error("Swept parameters must be bounded");

        Real prlb = bdcs.lowerBound[compIdx], prub = bdcs.upperBound[compIdx];
        Real lb, ub;

        if (args.count("range")) {
            auto rangeStr = args["range"].as<string>();
            vector<string> rangeComponents;
            boost::trim(rangeStr), boost::split(rangeComponents, rangeStr, boost::is_any_of(":"));
            if (rangeComponents.size() != 2) throw runtime_error("Invalid range; expected lower:upper");
            lb = stod(rangeComponents[0]), ub = stod(rangeComponents[1]);
        }
        else {
            Real rr = args["rangeRelative"].as<double>();
            Real prSize = prub - prlb;
            lb = params[compIdx] - rr * prSize, ub = params[compIdx] + rr * prSize;
        }

        if ((lb < prlb) || (ub > prub)) {
            std::cerr << "WARNING: Specified sweep range of " << lb << ":" << ub
                << " outside parameter range of " << prlb << ":" << prub << std::endl;
        }

        params[compIdx] = lb; // make dummy iterate below actually 0th iterate.
        cout << "it\tparam\tJFull\tgradp JFull";
        {
            auto &it = imanager->get(params.size(), params.data());
            for (const auto &etermptr : it.evaluatedObjectiveTerms())
                cout << "\t" << etermptr->name << "\tgradp " << etermptr->name;
        }
        cout << endl;

        const size_t nsamples = args["nsamples"].as<size_t>();
        for (size_t i = 0; i < nsamples; ++i) {
            if (args.count("singleIteration")) i = args["singleIteration"].as<size_t>();
            params[compIdx] = lb + ((nsamples == 1) ? 0.0 : (ub - lb) * (double(i) / (nsamples - 1)));
            auto &it = imanager->get(params.size(), params.data());
            cout << i << "\t" << params[compIdx] << "\t" << it.evaluate() << "\t" << it.gradp()[compIdx];
            for (const auto &etermptr : it.evaluatedObjectiveTerms())
                cout << "\t" << etermptr->value() << "\t" << etermptr->gradp[compIdx];
            cout << endl;

            if (args.count("output")) it.writeMeshAndFields(args["output"].as<string>() + "_" + std::to_string(i) + ".msh");
            if (args.count("singleIteration")) break;
        }

        // BENCHMARK_REPORT();
        return;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Run the optimizer
    ////////////////////////////////////////////////////////////////////////////
    string solver = args["solver"].as<string>(), output;
    if (args.count("output")) output = args["output"].as<string>();

    PO::OptimizerConfig oconfig;
    if (args.count("nIters")) oconfig.niters = args["nIters"].as<size_t>();
    oconfig.gd_step = args["step"].as<double>();

    if (solver == "lbfgs") oconfig.lbfgs_memory = 10;
    optimizers.at(solver)(params, bdcs, *imanager, oconfig, output);

    ////////////////////////////////////////////////////////////////////////////
    // Extract and process the result.
    ////////////////////////////////////////////////////////////////////////////
    std::vector<Real> result(params.domainSize());
    for (size_t i = 0; i < result.size(); ++i)
        result[i] = params[i];

    if (inflator.isParametric()) {
        cout << "Final p:";
        for (size_t i = 0; i < result.size(); ++i)
            cout << "\t" << result[i];
        cout << endl;
    }

    BENCHMARK_REPORT();
}

int main(int argc, const char *argv[])
{
    cout << setprecision(16);

    po::variables_map args = parseCmdLine(argc, argv);

    auto job = PO::parseJobFile(args["job"].as<string>());

    size_t deg = args["degree"].as<size_t>();
    if (auto job2D = dynamic_cast<PO::Job<2> *>(job)) {
        if (deg == 1) execute<2, 1>(args, job2D);
        if (deg == 2) execute<2, 2>(args, job2D);
    }
    else if (auto job3D = dynamic_cast<PO::Job<3> *>(job)) {
        if (deg == 1) execute<3, 1>(args, job3D);
        if (deg == 2) execute<3, 2>(args, job3D);
    }
    else throw std::runtime_error("Invalid job file.");

    return 0;
}
