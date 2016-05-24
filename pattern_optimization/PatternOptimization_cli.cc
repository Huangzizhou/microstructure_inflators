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
#include "LinearElasticity.hh"
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include "GlobalBenchmark.hh"

#include "Inflator.hh"
#include "MakeInflator.hh"

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <map>
#include <functional>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "optimizers/ceres.hh"
#include "optimizers/dlib.hh"
#include "optimizers/gradient_descent.hh"

#include "PatternOptimizationIterate.hh"

#include "PatternOptimizationJob.hh"
#include "BoundConstraints.hh"
#include "IterateFactory.hh"
#include "IterateManager.hh"

#include "OptimizerConfig.hh"
#include "PatternOptimizationConfig.hh"
#include "objective_terms/TensorFit.hh"
#include "objective_terms/ProximityRegularization.hh"

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;

using OptimizerMap =
    map<string, std::function<void(ScalarField<Real> &, const BoundConstraints &,
            IterateManagerBase &, const OptimizerConfig &, const string &)>>;
OptimizerMap optimizers = {
    {"levenberg_marquardt",  optimize_ceres_lm},
    {"dogleg",               optimize_ceres_dogleg},
    {"bfgs",                 optimize_dlib_bfgs},
    {"lbfgs",                optimize_dlib_bfgs},
    {"gradient_descent",     optimize_gd}
};

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
    visible_opts.add_options()("help",        "Produce this help message")
        ("inflator,i",   po::value<string>()->default_value("Isosurface"),       "inflator to use (defaults to Isosurface)")
        ("pattern,p",    po::value<string>(),                                    "Pattern wire mesh (.obj|wire)")
        ("material,m",   po::value<string>(),                                    "base material")
        ("degree,d",     po::value<size_t>()->default_value(2),                  "FEM Degree")
        ("output,o",     po::value<string>()->default_value(""),                 "output .js mesh + fields at each iteration")
        ("cell_size,c",  po::value<double>(),                                    "Inflation cell size (James' inflator only. Default: 5mm)")
        ("isotropicParameters,I",                                                "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                                    "Use vertex thickness instead of edge thickness (3D only)")
        ("subdivide,S",  po::value<size_t>()->default_value(0),                  "number of subdivisions to run for 3D inflator")
        ("sub_algorithm,A", po::value<string>()->default_value("simple"),        "subdivision algorithm for 3D inflator (simple or loop)")
        ("max_volume,v", po::value<double>(),                                    "maximum element volume parameter for wire inflator")
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs, levenberg_marquardt")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient descent step size")
        ("nIters,n",     po::value<size_t>(),                                    "number of iterations (infinite by default)")
        ("fullCellInflator",                                                     "use the full periodic inflator instead of the reflection-based one")
        ("ignoreShear",                                                          "Ignore the shear components in the isotropic tensor fitting")
        ("proximityRegularizationWeight", po::value<double>(),                   "Use a quadratic proximity regularization term with the specified weight.")
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

    size_t d = vm["degree"].as<size_t>();
    if (d < 1 || d > 2) {
        cout << "Error: FEM Degree must be 1 or 2" << endl;
        fail = true;
    }

    if (optimizers.count(vm["solver"].as<string>()) == 0) {
        cout << "Illegal solver specified" << endl;
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
void execute(const po::variables_map &args, const Job<_N> *job)
{
    auto infl_ptr = make_inflator<_N>(args["inflator"].as<string>(),
                                     filterInflatorOptions(args),
                                     job->parameterConstraints);

    Inflator<_N> &inflator = *infl_ptr;

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    cout << "target tensor: " << targetC << endl;

    SField params = job->validatedInitialParams(inflator);

    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;

    // Set up simulators' (base) material
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());

    using Iterate = Iterate<Simulator>;
    auto ifactory = make_iterate_factory<Iterate,
                                         ObjectiveTerms::IFConfigTensorFit<Simulator>,
                                         ObjectiveTerms::IFConfigProximityRegularization>(inflator);

    ifactory->ignoreShear = args.count("ignoreShear");
    if (ifactory->ignoreShear) cout << "Ignoring shear components" << endl;
    ifactory->ObjectiveTerms::IFConfigProximityRegularization::enabled = false;
    if (args.count("proximityRegularizationWeight")) {
        ifactory->ObjectiveTerms::IFConfigProximityRegularization::enabled    = true;
        ifactory->ObjectiveTerms::IFConfigProximityRegularization::initParams = job->initialParams;
        ifactory->ObjectiveTerms::IFConfigProximityRegularization::weight     = args["proximityRegularizationWeight"].as<double>();
    }

    auto imanager = make_iterate_manager(ifactory);
    ifactory->targetS = targetS;
    BoundConstraints bdcs(inflator, job->radiusBounds, job->translationBounds, job->blendingBounds,
                          job->varLowerBounds, job->varUpperBounds);

    string solver = args["solver"].as<string>(),
           output = args["output"].as<string>();
    OptimizerConfig oconfig;
    if (args.count("nIters")) oconfig.niters = args["nIters"].as<size_t>();
    oconfig.gd_step = args["step"].as<double>();

    if (solver == "lbfgs") oconfig.lbfgs_memory = 10;
    optimizers.at(solver)(params, bdcs, *imanager, oconfig, output);

    if (inflator.isParametric()) {
        cout << "Final p:";
        vector<Real> result(params.domainSize());
        for (size_t i = 0; i < result.size(); ++i) {
            result[i] = params[i];
            cout << "\t" << params[i];
        }
        cout << endl;
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
    else throw runtime_error("Invalid job file.");

    return 0;
}
