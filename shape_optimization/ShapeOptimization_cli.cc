////////////////////////////////////////////////////////////////////////////////
// ShapeOptimization_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      General shape optimization solver.
//      User should provide executable with wire, initial parameters, boundary
//      conditions and other options. Result should be the optimized shape (e.g.
//      the one with lower maximum stress for the given boudary conditions).
//      Strongly based on WCSOptimization_cli.cc
//
//      Example invocation:
//      ./ShapeOptimization_cli -p square.wire -b bondary_conditions.json
//              -m $MICRO_DIR/materials/B9Creator.material test_2D_job.opt
//
*/
//  Author:  Davi Colli Tozoni (dctozoni) davi.tozoni@nyu.edu
//  Company:  New York University
//  Created:  1/9/18
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <Materials.hh>
#include <GlobalBenchmark.hh>
#include <Future.hh>

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

#include <json.hpp>

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
#include <optimizers/gradient_descent.hh>
#include <optimizers/nlopt.hh>

#include <ShapeOptimizationIterate.hh>

#include <OptimizerConfig.hh>
#include <PatternOptimizationConfig.hh>
#include <objective_terms/ProximityRegularization.hh>
#include "StressObjectiveTerm.hh"

#include <constraints/Printability.hh>

#include <Parallelism.hh>

namespace po = boost::program_options;
namespace PO = PatternOptimization;
namespace SO = ShapeOptimization;
using json = nlohmann::json;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: ShapeOptimization_cli [options] job.opt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

using OptimizerMap =
map<string, std::function<void(ScalarField<Real> &, const PO::BoundConstraints &,
                               PO::IterateManagerBase &, const PO::OptimizerConfig &, const string &)>>;

OptimizerMap optimizers = {
        {"levenberg_marquardt",  optimize_ceres_lm},
        {"slsqp",                optimize_nlopt_slsqp},
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
            ("pattern,p",    po::value<string>(),                              "Pattern wire mesh (.obj|wire), or initial mesh for BoundaryPerturbationInflator")
            ("params",       po::value<string>(),                              "Initial params (overrides those specified in job file).")
            ;

    po::options_description simulationOptions;
    patternOptions.add_options()
            ("boundaryConditions,b", po::value<string>(),                    "boundary conditions")
            ;

    po::options_description meshingOptions;
    meshingOptions.add_options()
            ("meshingOptions,M",      po::value<string>(), "Meshing options configuration file")
            ;

    po::options_description optimizerOptions;
    optimizerOptions.add_options()
            ("nIters,n",     po::value<size_t>()->default_value(20),                 "number of iterations")
            ("step,s",       po::value<double>()->default_value(0.0001),             "gradient step size")
            ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, slsqp, levenberg_marquardt")
            ;

    po::options_description objectiveOptions;
    objectiveOptions.add_options()
            ("pnorm,P",      po::value<double>()->default_value(1.0),         "pnorm used in the Lp stress measure")
            ("usePthRoot,R",                                                  "Use the true Lp norm for global worst case stress measure (applying pth root)")
            ("stressWeight",    po::value<double>()->default_value(1.0),         "Weight for the Microscopic stress term of the objective")
            ("proximityRegularizationWeight", po::value<double>(),            "Use a quadratic proximity regularization term with the specified weight.")
            ("proximityRegularizationTarget", po::value<string>(),            "The target parameter values for the proximity regularization term (defaults to initial parameters.)")
            ;

    po::options_description constraintOptions;
    constraintOptions.add_options()
            ("PrintabilityConstraint", "Enforce self-supporting printability constraints as inequality constraints (for optimizers that support this)")
            ;

    po::options_description elasticityOptions;
    elasticityOptions.add_options()
            ("material,m",   po::value<string>(),                    "Base material")
            ("degree,d",     po::value<size_t>()->default_value(2),  "FEM Degree")
            ;

    po::options_description generalOptions;
    generalOptions.add_options()
            ("help,h",                                    "Produce this help message")
            ("output,o",             po::value<string>(), "Output mesh and fields at each iteration")
            ("dumpShapeDerivatives", po::value<string>(), "Dump shape derivative fields for Stress")
            ("numProcs",             po::value<size_t>(), "Number of threads to use for TBB parallelism (CGAL mesher, etc.)")
            ;

    po::options_description visibleOptions;
    visibleOptions.add(patternOptions).add(simulationOptions).add(meshingOptions)
            .add(optimizerOptions).add(objectiveOptions).add(constraintOptions)
            .add(elasticityOptions).add(generalOptions);

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
        if (vm["inflator"].as<string>().compare("HexaPillars")) {
            cout << "Error: must specify pattern mesh" << endl;
            fail = true;
        }
    }

    size_t d = vm["degree"].as<size_t>();
    if (d < 1 || d > 2) {
        cout << "Error: FEM Degree must be 1 or 2" << endl;
        fail = true;
    }

    if (optimizers.count(vm["solver"].as<string>()) == 0) {
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
void execute(const po::variables_map &args, PO::Job<_N> *job)
{
    auto infl_ptr = make_inflator<_N>(args["inflator"].as<string>(),
                                      filterInflatorOptions(args),
                                      job->parameterConstraints);
    auto &inflator = *infl_ptr;

    string bcondsPath = args["boundaryConditions"].as<string>();

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    // Parse parameters
    auto parseParams = [](string pstring) -> vector<Real> {
        boost::trim(pstring);
        vector<string> tokens;
        boost::split(tokens, pstring, boost::is_any_of("\t "),
                     boost::token_compress_on);
        vector<Real> pvals;
        for (string &s : tokens) pvals.push_back(std::stod(s));
        return pvals;
    };

    // If requested, override the initial parameters set in the job file
    if (args.count("params"))
        job->initialParams = parseParams(args["params"].as<string>());

    SField params = job->validatedInitialParams(inflator);

    PO::BoundConstraints bdcs(inflator, job->radiusBounds, job->translationBounds, job->blendingBounds, job->metaBounds,
                              job->custom1Bounds, job->custom2Bounds, job->custom3Bounds, job->custom4Bounds,
                              job->custom5Bounds, job->custom6Bounds, job->custom7Bounds, job->custom8Bounds,
                              job->varLowerBounds, job->varUpperBounds);

    using StressTermConfig    = PO::ObjectiveTerms::IFConfigMicroscopicStress<Simulator>;
    using PRegTermConfig      = PO::ObjectiveTerms::IFConfigProximityRegularization;
    using PConstraintConfig   = PO::   Constraints::IFConfigPrintability<Simulator>;

    auto ifactory = PO::make_iterate_factory<SO::Iterate<Simulator>,
            StressTermConfig,
            PRegTermConfig,
            PConstraintConfig>(inflator, bdcs);

    ////////////////////////////////////////////////////////////////////////////
    // Configure the objective terms
    ////////////////////////////////////////////////////////////////////////////
    ifactory->StressTermConfig   ::enabled = args["StressWeight"].as<double>() != 0;
    ifactory->PRegTermConfig     ::enabled = args.count("proximityRegularizationWeight");
    ifactory->PConstraintConfig  ::enabled = false;

    // Configure Stress Objective
    // By default, an "Lp norm" objective is really the p^th power of the Lp norm.
    // To use the true "Lp norm", globalObjectiveRoot must be set to
    // 2.0 * globalObjectivePNorm (since pointwise Stress is already squared (e.g. Frobenius) norm)
    ifactory->StressTermConfig::weight = args["StressWeight"].as<double>();
    Real pnorm = args["pnorm"].as<double>();
    ifactory->StressTermConfig::globalObjectivePNorm = pnorm;
    ifactory->StressTermConfig::globalObjectiveRoot  = args.count("usePthRoot") ? 2.0 * pnorm : 1.0;

    if (args.count("PrintabilityConstraint")) {
        ifactory->PConstraintConfig::enabled = true;
    }

    if (args.count("proximityRegularizationWeight")) {
        ifactory->PRegTermConfig::enabled = true;
        ifactory->PRegTermConfig::weight = args["proximityRegularizationWeight"].as<double>();
        ifactory->PRegTermConfig::targetParams = job->initialParams;
        if (args.count("proximityRegularizationTarget")) {
            ifactory->PRegTermConfig::targetParams = parseParams(args["proximityRegularizationTarget"].as<string>());
            if (ifactory->PRegTermConfig::targetParams.size() != job->initialParams.size())
                throw runtime_error("Invalid proximity regularization target parameter count");
        }
        if (args.count("proximityRegularizationZeroTarget")) {
            vector<Real> zeroTargetParams(job->numParams(), 0.0);
            ifactory->PRegTermConfig::IFConfigProximityRegularization::targetParams = zeroTargetParams;
        }
    }

    auto imanager = PO::make_iterate_manager(std::move(ifactory));

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

    json output_json;
    if (inflator.isParametric()) {
        output_json["final_p"] = json::array();
        cout << "Final p:";
        for (size_t i = 0; i < result.size(); ++i) {
            cout << "\t" << result[i];
            output_json["final_p"].push_back(result[i]);
        }
        cout << endl;
    }

    BENCHMARK_REPORT_NO_MESSAGES();
}

int main(int argc, const char *argv[]) {
    po::variables_map args = parseCmdLine(argc, argv);

#if HAS_TBB
    size_t np = tbb::task_scheduler_init::default_num_threads();
    if (args.count("numProcs")) {
        size_t manualNP = args["numProcs"].as<size_t>();
        if (manualNP > np)
            std::cerr << "WARNING: specifying more than the default number of TBB threads." << std::endl;
        np = manualNP;
    }
    tbb::task_scheduler_init init(np);
#else
    if (args.count("numProcs"))
        std::cerr << "WARNING: parallelism disabled; numProcs argument ignored." << std::endl;
#endif

    cout << setprecision(16);
    auto job = PO::parseJobFile(args["job"].as<string>());

    size_t deg = args["degree"].as<size_t>();
    if (auto job2D = dynamic_cast<PO::Job<2> *>(job.get())) {
        if (deg == 1) execute<2, 1>(args, job2D);
        if (deg == 2) execute<2, 2>(args, job2D);
    }
    else if (auto job3D = dynamic_cast<PO::Job<3> *>(job.get())) {
        if (deg == 1) execute<3, 1>(args, job3D);
        if (deg == 2) execute<3, 2>(args, job3D);
    }
    else throw std::runtime_error("Invalid job file.");

    return 0;
}


