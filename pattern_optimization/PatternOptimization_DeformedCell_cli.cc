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
#include <boost/algorithm/string.hpp> // requiured for parsing jacobian


#include "PatternOptimization.hh"
#include "PatternOptimizationJob.hh"

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;
using namespace PeriodicHomogenization;


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
        ("jacobian,j",   po::value<string>()->default_value("1.0 0.0 0.0 1.0"),  "linear deformation")
        ("final_mesh,f", po::value<string>(), "output .msh file name prefix")
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

    if (vm.count("final_mesh") == 0) {
        cout << "Error: must specify output mesh name prefix" << endl;
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
void execute_defCell(const po::variables_map args,
                     const vector<MeshIO::IOVertex> inVertices,
                     const vector<MeshIO::IOElement> inElements,
                     const Job<_N> *job)
{
    auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());
    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;
    typedef typename Simulator::VField VField;
    Simulator sim(inElements, inVertices);
    auto &mesh = sim.mesh();

    // Parse jacobian.
    vector<string> jacobianComponents;
    string jacobianString = args["jacobian"].as<string>();
    boost::trim(jacobianString);
    boost::split(jacobianComponents, jacobianString, boost::is_any_of("\t "),
                 boost::token_compress_on);
    if (jacobianComponents.size() != _N * _N)
        throw runtime_error("Invalid deformation jacobian");
    Eigen::Matrix<Real, _N, _N> jacobian;
    for (size_t i = 0; i < _N; ++i) {
        for (size_t j = 0; j < _N; ++j) {
            jacobian(i, j) = stod(jacobianComponents[_N * i + j]);
        }
    }

    auto bbox = mesh.boundingBox();
    VectorND<_N> center = 0.5 * (bbox.minCorner + bbox.maxCorner);
    vector<MeshIO::IOVertex> deformedVertices;

    for (size_t vi = 0; vi < mesh.numVertices(); ++vi) {
        VectorND<_N> p = mesh.vertex(vi).node()->p;
        deformedVertices.emplace_back((jacobian * (p - center) + center).eval());
    }

    //Real deformedCellVolume = bbox.volume() * jacobian.determinant();

    vector<VField> w_ij;
    // Morteza's transformation formulas
    mat.setTensor(mat.getTensor().transform(jacobian.inverse()));
    solveCellProblems(w_ij, sim);
    auto EhDefo = homogenizedElasticityTensor(w_ij, sim).transform(jacobian);
    cout << setprecision(16);
    cout << "Elasticity tensor:" << endl;
    cout << EhDefo << endl << endl;
    cout << "and its anisotropy is: " << EhDefo.anisotropy() << endl;
    cout << "Homogenized Moduli: ";
    EhDefo.printOrthotropic(cout);
    MeshIO::save(args["final_mesh"].as<string>()+"_ref.msh", inVertices, inElements);
    MeshIO::save(args["final_mesh"].as<string>()+"_def.msh", deformedVertices, inElements);


	// printing out the target C
    auto targetC = job->targetMaterial.getTensor();

	cout << "target Moduli was: " << endl;
	cout << targetC << endl;
	cout << "and its anisotropy is: " << targetC.anisotropy() << endl;
}


template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args, const Job<_N> *job)
{
    shared_ptr<ConstrainedInflator<_N>> inflator_ptr;
    if (_N == 2) {
        inflator_ptr = make_shared<ConstrainedInflator<_N>>(
                job->parameterConstraints,
                args["pattern"].as<string>(),
                3); // MHS JUL14, 2015: the last parameter is a symmetryMode (see the corresponding constructor in Inflator.hh for more details)
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


    std::cout << "numer of inflator parametrs is " << inflator.numParameters() << endl;

    /////*** things needed for deformed cell optimization ***/////
    // parse jacobian
    vector<string> jacobianComponents;
    string jacobianString = args["jacobian"].as<string>();
    boost::trim(jacobianString);
    boost::split(jacobianComponents, jacobianString, boost::is_any_of("\t "),
                 boost::token_compress_on);
    if (jacobianComponents.size() != _N * _N)
        throw runtime_error("Invalid deformation jacobian");
    Eigen::Matrix<Real, _N, _N> jacobian;
    for (size_t i = 0; i < _N; ++i) {
        for (size_t j = 0; j < _N; ++j) {
            jacobian(i, j) = stod(jacobianComponents[_N * i + j]);
        }
    }


    // use the deformed target C
    auto targetC = job->targetMaterial.getTensor().transform(jacobian.inverse());
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

    // Morteza's transformation formulas
    mat.setTensor(mat.getTensor().transform(jacobian.inverse()));

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

    cout << "writing down the undeformed and deformed final meshes." << endl;
    execute_defCell<_N,_FEMDegree>(args, inflator.vertices(), inflator.elements(), job);

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