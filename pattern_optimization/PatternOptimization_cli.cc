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
    visible_opts.add_options()("help",        "Produce this help message")
        ("pattern,p",    po::value<string>(), "Pattern wire mesh (.obj|wire)")
        ("material,m",   po::value<string>(), "base material")
        ("degree,d",     po::value<size_t>()->default_value(2),                  "FEM Degree")
        ("output,o",     po::value<string>()->default_value(""),                 "output .js mesh + fields at each iteration")
        ("dofOut",       po::value<string>()->default_value(""),                 "output pattern dofs in James' format at each iteration (3D Only)")
        ("cell_size,c",  po::value<double>()->default_value(5.0),                "Inflation cell size (3D only)")
        ("isotropicParameters,I",                                                "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                                    "Use vertex thickness instead of edge thickness (3D only)")
        ("subdivide,S",  po::value<size_t>()->default_value(0),                  "number of subdivisions to run for 3D inflator")
        ("sub_algorithm,A", po::value<string>()->default_value("simple"),        "subdivision algorithm for 3D inflator (simple or loop)")
        ("max_volume,v", po::value<double>(),                                    "maximum element volume parameter for wire inflator")
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs, levenberg_marquardt")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient step size")
        ("nIters,n",     po::value<size_t>()->default_value(20),                 "number of iterations")
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

    set<string> solvers = {"gradient_descent", "bfgs", "lbfgs", "levenberg_marquardt", "lm_bd_penalty", "levmar"};
    if (solvers.count(vm["solver"].as<string>()) == 0) {
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
    shared_ptr<ConstrainedInflator<_N>> inflator_ptr;
    if (_N == 2) {
        inflator_ptr = make_shared<ConstrainedInflator<_N>>(
                job->parameterConstraints,
                args["pattern"].as<string>());
    }
    else {
        inflator_ptr = make_shared<ConstrainedInflator<_N>>(
                job->parameterConstraints,
                args["pattern"].as<string>(),
                args["cell_size"].as<double>(),
                0.5 * sqrt(2),
                args.count("isotropicParameters"),
                args.count("vertexThickness"));
        inflator_ptr->configureSubdivision(args["sub_algorithm"].as<string>(),
                                           args["subdivide"].as<size_t>());
    }

	ConstrainedInflator<_N> &inflator = *inflator_ptr;
    if (args.count("max_volume"))
        inflator.setMaxElementVolume(args["max_volume"].as<double>());

    auto targetC = job->targetMaterial.getTensor();
    ETensor<_N> targetS = targetC.inverse();

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    cout << "target tensor: " << targetC << endl;

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

    string dofOut = args["dofOut"].as<string>();
    if (dofOut != "")
        inflator.setDoFOutputPrefix(dofOut);

    SField params(job->initialParams);
    for (const auto &boundEntry : job->varLowerBounds) {
        if (boundEntry.first > params.domainSize())
            cerr << "WARNING: bound on nonexistent variable" << endl;
    }

    for (size_t p = 0; p < params.domainSize(); ++p) {
        if (job->varLowerBounds.count(p)) {
             if ((params[p] < job->varLowerBounds.at(p)) ||
                 (params[p] > job->varUpperBounds.at(p))) {
                throw std::runtime_error("Initial point infeasible");
             }
        }
    }

    Optimizer<Simulator> optimizer(inflator, job->radiusBounds, job->translationBounds,
                                   job->varLowerBounds, job->varUpperBounds);
    string solver = args["solver"].as<string>(),
           output = args["output"].as<string>();
    size_t niters = args["nIters"].as<size_t>();
    if (solver == "levenberg_marquardt")
        optimizer.optimize_lm(params, targetS, output);
    else if (solver == "lm_bd_penalty")
        optimizer.optimize_lm_bound_penalty(params, targetS, output);
    else if (solver == "levmar")
        optimizer.optimize_levmar(params, targetS, niters, output);
    else if (solver == "gradient_descent")
        optimizer.optimize_gd(params, targetS, niters,
                          args["step"].as<double>(), output);
    else if (solver == "bfgs")
        optimizer.optimize_bfgs(params, targetS, niters, output);
    else if (solver == "lbfgs")
        optimizer.optimize_bfgs(params, targetS, niters, output, 10);

    std::cout << "Final p:";
    std::vector<Real> result(params.domainSize());
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = params[i];
        cout << "\t" << params[i];
    }
    if (dofOut != "") inflator.writePatternDoFs(dofOut + ".final.dof", result);
    cout << endl;


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
    else throw std::runtime_error("Invalid job file.");

    return 0;
}
