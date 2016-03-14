////////////////////////////////////////////////////////////////////////////////
// PatternOptimization_DeformedCell_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Evolves the microstructure parameters to bring the structure's
//      homogenized elasticity tensor of the deformed closer to a target tensor.
*/ 
//  Author:  Morteza H Siboni, hakimi1364@gmail.com
//  Note  :  Base on Julinan's PatternOptmization_cli.cc
//  Company:  New York University
//  Created:  2015 
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include "LinearElasticity.hh"
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include "GlobalBenchmark.hh"
#include "EdgeMeshType.h"
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
    cout << "Usage: PatternOptimization_DeformedCell_cli [options] job.opt" << endl;
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
        ("jacobians,j",  po::value<string>(), "a .txt file includeing the jacobinas 'J11 J12 J21 J22'")
        ("sym",          po::value<int>()->default_value(3), "symmetry mode, use -1 to fall back to Luigi's symmetry mode")
        ("regularization,r",          po::value<double>()->default_value(0.0), "regularization weight")
        ("degree,d",     po::value<size_t>()->default_value(2),                  "FEM Degree")
        ("output,o",     po::value<string>()->default_value(""),                 "output .js mesh + fields at each iteration")
        ("dofOut",       po::value<string>()->default_value(""),                 "output pattern dofs in James' format at each iteration (3D Only)")
        ("cell_size,c",  po::value<double>()->default_value(5.0),                "Inflation cell size (3D only)")
        ("isotropicParameters,I",                                                "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                                    "Use vertex thickness instead of edge thickness (3D only)")
        ("subdivide,S",  po::value<size_t>()->default_value(0),                  "number of subdivisions to run for 3D inflator")
        ("sub_algorithm,A", po::value<string>()->default_value("simple"),        "subdivision algorithm for 3D inflator (simple or loop)")
        ("max_volume,v", po::value<double>(),                                    "maximum element volume parameter for wire inflator")
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs, levenberg_marquardt, levenberg_marquardt_reg")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient step size")
        ("nIters,n",     po::value<size_t>()->default_value(20),                 "number of iterations")
        ("fullCellInflator",                                                     "use the full periodic inflator instead of the reflection-based one")
        ("patternOut",   po::value<string>()->default_value(""),                                    "filename.txt that includes the optimized pattern parameters on a square cell 'p1 p2 ...'")
        ("stiffnessOut", po::value<string>()->default_value(""),                                    "filename.txt that includes the optimized stiffness in flatened format 'C11 ...'")
        ("costOut",      po::value<string>()->default_value(""),                                    "filename.txt that includes the intial/final/reltive change of the cost function  'initial final costChange'")
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

    set<string> solvers = {"gradient_descent", "bfgs", "lbfgs", "levenberg_marquardt", "levenberg_marquardt_reg", "lm_bd_penalty", "levmar"};
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


	// targetC - C
	cout << "printing C-Diff:" << endl;
	ETensor<_N> finC = EhDefo;
	ETensor<_N> tarC = targetC;
	ETensor<_N> diffC = finC - tarC;

	cout << diffC << endl;
	cout << "its norm is: " << diffC.quadrupleContract(diffC) / 2.0 << endl;
	

	// targetS - S
	cout << "printing S-Diff:" << endl;
	ETensor<_N> finS = EhDefo.inverse();
	ETensor<_N> tarS = targetC.inverse();
	ETensor<_N> diffS = finS - tarS;

	cout << diffS << endl;
	cout << "its norm is: " << diffS.quadrupleContract(diffS) / 2.0 << endl;
	
}

template<size_t _N>
void readDeformations(const po::variables_map &args, vector<Eigen::Matrix<Real, _N, _N>> &defs)
{
	defs.clear();
	ifstream in;
	in.open(args["jacobians"].as<string>());

	Real J11, J12, J21, J22;
	while(in >> J11 >> J12 >> J21 >> J22)
	{
		Eigen::Matrix<Real, _N, _N> jacobian;
		jacobian << J11, J12,
				    J21, J22;
		defs.push_back(jacobian);
	}

}

// new execute TODO edit it 
template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args, 
			 const Job<_N> *job, 
			 const Eigen::Matrix<Real, _N, _N> deformation, 
			 vector<MeshIO::IOVertex>  & refMeshVertices,
             vector<MeshIO::IOElement> & refMeshElements,
			 vector<Real>            & outParams,
			 vector<Real>            & stiffness,
			 vector<Real>            & costs) 
{

    shared_ptr<ConstrainedInflator<_N>> inflator_ptr;
    if (_N == 2) {

		if (args["sym"].as<int>() >= -1 && args["sym"].as<int>() < 8)
			inflator_ptr = make_shared<ConstrainedInflator<_N>>(
					job->parameterConstraints,
					args["pattern"].as<string>(),
					args["sym"].as<int>()); // this  parameter is a symmetryMode (see the corresponding constructor in Inflator.hh for more details)
		else
			throw("symmetry mode must be in [-1..7]");

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
        inflator_ptr->setReflectiveInflator(args.count("fullCellInflator") == 0);
    }

	ConstrainedInflator<_N> &inflator = *inflator_ptr;
    if (args.count("max_volume"))
        inflator.setMaxElementVolume(args["max_volume"].as<double>());


    // use the deformed target C
    auto targetC = job->targetMaterial.getTensor().transform(deformation.inverse());

    ETensor<_N> targetS = targetC.inverse();

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
    mat.setTensor(mat.getTensor().transform(deformation.inverse()));

    string dofOut = args["dofOut"].as<string>();
    if (dofOut != "")
        inflator.setDoFOutputPrefix(dofOut);


	std::cout << "this is the target C inputed to the optimizer :" << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
	targetC.printOrthotropic(std::cout);
	std::cout << "-----------------------------------------------" << std::endl;


	SField params(job->initialParams);
	SField initialParams(job->initialParams);


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

    Real initialCost, finalCost;
    ETensor<_N> transformedOptimizedC;

    Optimizer<Simulator> optimizer(inflator, job->radiusBounds, job->translationBounds,
                                   job->varLowerBounds, job->varUpperBounds);


    string output = args["output"].as<string>();

	string solver = args["solver"].as<string>();
    size_t niters = args["nIters"].as<size_t>();

    if (solver == "levenberg_marquardt")
        optimizer.optimize_lm(params, targetS, output);
    else if (solver == "levenberg_marquardt_reg")
        optimizer.optimize_lm_regularized(params, initialParams, args["regularization"].as<double>(), targetS, output, initialCost, finalCost, transformedOptimizedC);
    else if (solver == "lm_bd_penalty")
        optimizer.optimize_lm_bound_penalty(params, targetS, output);
    else if (solver == "levmar")
        optimizer.optimize_levmar(params, targetS, niters, output);
    else if (solver == "gradient_descent")
        optimizer.optimize_gd(params, targetS, niters,
                          args["step"].as<double>(), output);
    else if (solver == "bfgs")
        optimizer.optimize_bfgs(params, targetS, niters, output);
	/* else if (solver == "bfgs_reg") */
        /* optimizer.optimize_bfgs_regularized(params, initialParams, args["regularization"].as<double>(), targetS, output1, output2,  niters, 10); */
    else if (solver == "lbfgs")
        optimizer.optimize_bfgs(params, targetS, niters, output, 10);

    //std::cout << "Final p:";
    std::vector<Real> result(params.domainSize());
    outParams.clear();
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = params[i];
        //cout << "\t" << params[i];
        outParams.push_back(params[i]);
    }


	costs.push_back(initialCost);
	costs.push_back(finalCost);
	costs.push_back(std::abs(finalCost - initialCost));

	ETensor<_N> C = transformedOptimizedC.transform(deformation);
	std::cout << "the anisotopy of final C [computed outside the optimizer after applying the transformations] is  " << C.anisotropy() << std::endl;

	Real Ex, Ey, nuYX, muXY;
	C.getOrthotropic2D(Ex, Ey, nuYX, muXY);


	std::cout << "the orthotropic components  of final C are " << Ex << "\t"
															   << Ey << "\t" 
															   << nuYX << "\t" 
															   << muXY << std::endl;

	C.getUpperRight2D(stiffness);

    refMeshVertices = inflator.vertices();
    refMeshElements = inflator.elements();



    if (dofOut != "") inflator.writePatternDoFs(dofOut + ".final.dof", result);

    cout << endl;
}


template<size_t _N, size_t _FEMDegree>
void runPatternOptimization(const po::variables_map &args, 
							const Job<_N> *job, 
							const vector<Eigen::Matrix<Real, _N, _N>> defsTable, 
							vector<vector<Real>>           & paramsTable,
							vector<vector<Real>>           & stiffnessTable,
							vector<vector<Real>>           & costsTable)
{
	vector<MeshIO::IOVertex>   refMeshVertices;
	vector<MeshIO::IOElement>  refMeshElements;

	vector<Real> currentParams, currentStiffness, currentCosts;

	for (auto def = defsTable.begin(); def != defsTable.end(); ++def){
		execute<_N, _FEMDegree>(args, job, *def, refMeshVertices, refMeshElements, currentParams, currentStiffness, currentCosts);
		paramsTable.push_back(currentParams);
		stiffnessTable.push_back(currentStiffness);
		costsTable.push_back(currentCosts);
		
		currentParams.clear();
		currentStiffness.clear();
		currentCosts.clear();
	}
}

// helper functions to generate the output files/meshes
template<typename _Real>
void tableToFile(const vector<vector<_Real>> & table,
				 const std::string fileName)
{
	ofstream out;
	out.open(fileName);

	for (auto it = table.begin(); it != table.end(); ++it)
	{
		vector<_Real> row = *it;
		for (size_t i = 0; i < row.size() - 1; ++i)
		{
			out << showpos << scientific;
			out << row[i] << "\t";
		}
		out << row.back() << "\n";
	}
	out.close();
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


	vector<Eigen::Matrix<Real, 2, 2>> defs;
	readDeformations<2>(args, defs);

	vector<vector<Real>>	paramsTable;
	vector<vector<Real>>	stiffnessTable;
	vector<vector<Real>>	costsTable;

    size_t deg = args["degree"].as<size_t>();
	auto job2D = dynamic_cast<Job<2> *>(job);
	if (deg == 1) runPatternOptimization<2, 1>(args, job2D, defs, paramsTable, stiffnessTable, costsTable);
	if (deg == 2) runPatternOptimization<2, 2>(args, job2D, defs, paramsTable, stiffnessTable, costsTable);

	string patternOut = args["patternOut"].as<string>();
	if (patternOut != "")
		tableToFile<Real>(paramsTable, patternOut);

	string stiffnessOut = args["stiffnessOut"].as<string>();
	if (stiffnessOut != "")
		tableToFile<Real>(stiffnessTable, stiffnessOut);

	string costOut = args["costOut"].as<string>();
	if (costOut != "")
		tableToFile<Real>(costsTable, costOut);

    return 0;
}
