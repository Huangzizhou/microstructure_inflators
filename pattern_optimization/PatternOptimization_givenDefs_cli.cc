////////////////////////////////////////////////////////////////////////////////
//// PatternOptimization_quadMesh_cli.cc
//////////////////////////////////////////////////////////////////////////////////
///*! @file
////      This reads in a (2D) quad mesh, finds the equivalent parallelograms, and
////	  runs pattern optimization on them 
//*/ 
////  Author:  Morteza H Siboni (mhs), m.hakimi.siboni.@gmail.com
////  Company:  New York University
////  Created:  08/10/2014
//////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include "LinearElasticity.hh"
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include "GlobalBenchmark.hh"
// I thought these are needed to call static functions from WireMesh2D.h, but they are not!
//#include "EdgeMeshUtils.h"
//#include "WireMesh2D.h"

#include "PolyMeshType.h"
#include "PolyMeshUtils.h"
#include "WireMeshEmbedding.h"

#include "EdgeMeshType.h"
#include "Inflator.hh"

#include <vector>
#include <queue>
#include <cstdlib>
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
namespace fs = boost::filesystem;
using namespace fs;
using namespace std;
using namespace PatternOptimization;
using namespace PeriodicHomogenization;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: PatternOptimization_quadMesh_cli [options] job.opt" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}


po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("job", po::value<string>(), "job configuration file");

    po::positional_options_description p;
    p.add("job", 1);

    po::options_description visible_opts;
   	visible_opts.add_options()("help",			"Produce this help message")
        ("pattern,p",					po::value<string>(), "Pattern wire mesh (.obj|wire)")
        ("material,m",					po::value<string>(), "base material")
        ("quad,q",						po::value<string>(), "quad mesh (.obj|msh)")
//        ("out_folder,f",				po::value<string>(), "folder name to store the results")
        ("regularization,r",			po::value<double>()->default_value(0.0), "regularization weight")
//        ("trust,t",		 po::value<double>()->default_value(1e4), "trust region initial radius")
        ("degree,d",     po::value<size_t>()->default_value(2),                  "FEM Degree")
        ("output,o",     po::value<string>()->default_value(""),                 "output .js mesh + fields at each iteration")
        ("dofOut",       po::value<string>()->default_value(""),                 "output pattern dofs in James' format at each iteration (3D Only)")
        ("cell_size,c",  po::value<double>()->default_value(5.0),                "Inflation cell size (3D only)")
        ("isotropicParameters,I",                                                "Use isotropic DoFs (3D only)")
        ("vertexThickness,V",                                                    "Use vertex thickness instead of edge thickness (3D only)")
        ("subdivide,S",  po::value<size_t>()->default_value(0),                  "number of subdivisions to run for 3D inflator")
        ("sub_algorithm,A", po::value<string>()->default_value("simple"),        "subdivision algorithm for 3D inflator (simple or loop)")
        ("max_volume,v", po::value<double>(),                                    "maximum element volume parameter for wire inflator")
        ("solver",       po::value<string>()->default_value("gradient_descent"), "solver to use: none, gradient_descent, bfgs, lbfgs, levenberg_marquardt, dogleg")
        ("sym",          po::value<int>()->default_value(3), 				 "symmetry mode")
        ("step,s",       po::value<double>()->default_value(0.0001),             "gradient step size")
        ("nIters,n",     po::value<size_t>()->default_value(20),                 "number of iterations")
        ("fullCellInflator",                                                     "use the full periodic inflator instead of the reflection-based one")
        ;
	
	po::options_description cli_opts;
	cli_opts.add(visible_opts).add(hidden_opts);

	po::variables_map vm;

	try {
		po::store(po::command_line_parser(argc, argv).
				  options(cli_opts).positional(p).run(), vm);
		po::notify(vm);
	}
	catch (std::exception &e){
		cout << "Error: " << e.what() << endl << endl;
		usage(1, visible_opts);
	}

	bool fail = false;
	if (vm.count("job") == 0) {
		cout << "Error: must specify input job.opt file!" << endl;
		fail = true;
	}

	if (vm.count("pattern") == 0) {
		cout << "Error: must specify the edge mesh!" << endl;
		fail = true;
	}
		
	if (vm.count("material") == 0) {
		cout << "Error: must specify the material!" << endl;
		fail = true;
	}

	if (vm.count("quad") == 0) {
		cout << "Error: must specify the quad mesh!" << endl;
		fail = true;
	}
	
	/* if (vm.count("out_folder") == 0) { */
	/* 	cout << "Error: must specify a name for the output folder!" << endl; */
	/* 	fail = true; */
	/* } */

	if (fail || vm.count("help"))
		usage(fail, visible_opts);

	return vm;
}


template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;

// forward decleration of some of the helper functions
void tableToFile(const vector<ETensor<2>> & table, path savePath, const std::string fileName);
void tableToFile(const vector<vector<double>> & table, path savePath, const std::string fileName);
void tableToFile(const vector<Eigen::Matrix2d> & table, path savePath, const std::string fileName);

// read the quad mesh:
void readQuadMesh(const string & quadMeshPath, PolyMesh & pmesh)
{
	typedef WireMeshEmbedding<EMesh, PolyMesh>				WireEmbedding;
	typedef typename WireEmbedding::QuadParametrization 	QuadParametrization;
	typedef PolyMeshUtils<PolyMesh> PMU;
	bool ok = false;
	ok = PMU::importFromOBJ(quadMeshPath, pmesh);
	if (ok)
		WireEmbedding::preprocessQuadMesh(pmesh);

	// uncomment the following if you want to ignore Luigi's optimal parametrization
	/* for(auto fc = pmesh.face.begin(); fc != pmesh.face.end(); ++fc) */
	/* { */
	/* 	QuadParametrization & qpar = WireEmbedding::getQuadParametrizationHandle(pmesh)[fc]; */
	/* 	qpar.index0 = 0; */
	/* } */

	int faceCounter = 0;
	for(auto fc = pmesh.face.begin(); fc != pmesh.face.end(); ++fc)
	{
		cout << "points in face " << faceCounter << " are: " << "("  << fc->cP(0)[0] << "," << fc->cP(0)[1] << ") " <<
			                                                    "("  << fc->cP(1)[0] << "," << fc->cP(1)[1] << ") " <<
			                                                    "("  << fc->cP(2)[0] << "," << fc->cP(2)[1] << ") " <<
			                                                    "("  << fc->cP(3)[0] << "," << fc->cP(3)[1] << ")"  << endl; 
		
		QuadParametrization qpar = WireEmbedding::getQuadParametrizationHandle(pmesh)[fc];
		char index0 = qpar.index0;
		char index1 = (index0 + 1) % 4;
		char index2 = (index0 + 2) % 4;
		char index3 = (index0 + 3) % 4;
	
		cout << "points in face " << faceCounter << " are: " << "("  << fc->cP(index0)[0] << "," << fc->cP(index0)[1] << ") " <<
			                                                    "("  << fc->cP(index1)[0] << "," << fc->cP(index1)[1] << ") " <<
			                                                    "("  << fc->cP(index2)[0] << "," << fc->cP(index2)[1] << ") " <<
			                                                    "("  << fc->cP(index3)[0] << "," << fc->cP(index3)[1] << ")"  << endl; 

		++faceCounter;
	}
}

// get the deformation (F or U) for each quad in the quadMesh 
void getDeformations(PolyMesh & pmesh, 
		             vector<Eigen::Matrix2d> & deformations, 
		             vector<double> & angles,
		             const string type = "jacobian")
{
	typedef WireMeshEmbedding<EMesh, PolyMesh>				WireEmbedding;

	if (type == "jacobian")
		WireEmbedding::getJacobians(pmesh, deformations, angles);
	else if (type == "stretch")
		WireEmbedding::getStretches(pmesh, deformations, angles);
	else
		;
}


template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;

template<size_t _N>
using ETensor = ElasticityTensor<Real, _N>;
typedef ScalarField<Real> SField;


template<size_t _N, size_t _FEMDegree>
void execute(const po::variables_map &args, 
			 const Job<_N> *job, 
			 const Eigen::Matrix2d deformation, 
			 const int quadID,
			 vector<MeshIO::IOVertex>  & refMeshVertices,
             vector<MeshIO::IOElement> & refMeshElements,
			 vector<double>            & outParams) 
{
    shared_ptr<ConstrainedInflator<_N>> inflator_ptr;
    if (_N == 2) {
        inflator_ptr = make_shared<ConstrainedInflator<_N>>(
                job->parameterConstraints,
                args["pattern"].as<string>(),
                args["sym"].as<int>()); // this  parameter is a symmetryMode (see the corresponding constructor in Inflator.hh for more details)
        // this uses the original Luigi's symmetry stuff
        /* inflator_ptr = make_shared<ConstrainedInflator<_N>>( */
        /*         job->parameterConstraints, */
        /*         args["pattern"].as<string>()); // */
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


    std::cout << "numer of inflator parametrs is " << inflator.numParameters() << endl;


    // use the deformed target C
    auto targetC = job->targetMaterial.getTensor().transform(deformation.inverse());

    ETensor<_N> targetS = targetC.inverse();

    cout << "Target moduli:\t";
    targetC.printOrthotropic(cout);
    cout << endl;

    cout << "target tensor: " << endl<< targetC << endl;

  
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

	cout << "base material stiffness" << endl << mat.getTensor() << endl;

	/* char myC = 'c'; */
    /* while (myC != 'n') */
    	/* myC = getchar(); */

    string dofOut = args["dofOut"].as<string>();
    if (dofOut != "")
        inflator.setDoFOutputPrefix(dofOut);




	SField params(job->initialParams);
	SField initialParams(job->initialParams);


   	
	// cout << endl << params << endl;


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


    auto savePath = current_path() / args["output"].as<string>() / "inter_meshes";
    if (!exists(savePath) ||  !is_directory(savePath))
    	create_directory(savePath);
	
    string output1 = savePath.string();

	string quadFolder = "quad_" + to_string(quadID);
	savePath = current_path() / args["output"].as<string>() / "inter_meshes" / quadFolder;
    if (!exists(savePath) ||  !is_directory(savePath))
    	create_directory(savePath);

    string output2 = savePath.string();

    string output  = output1;

	string solver = args["solver"].as<string>();
    size_t niters = args["nIters"].as<size_t>();
    
    if (solver == "levenberg_marquardt")
        optimizer.optimize_lm(params, targetS, output);
    else if (solver == "levenberg_marquardt_reg")
        optimizer.optimize_lm_regularized(params, initialParams, args["regularization"].as<double>(), targetS, output1, output2);
    else if (solver == "lm_bd_penalty")
        optimizer.optimize_lm_bound_penalty(params, targetS, output);
    else if (solver == "levmar")
        optimizer.optimize_levmar(params, targetS, niters, output);
    else if (solver == "gradient_descent")
        optimizer.optimize_gd(params, targetS, niters,
                          args["step"].as<double>(), output);
    else if (solver == "bfgs")
        optimizer.optimize_bfgs(params, targetS, niters, output);
	else if (solver == "bfgs_reg")
        optimizer.optimize_bfgs_regularized(params, initialParams, args["regularization"].as<double>(), targetS, output1, output2,  niters, 10);
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

    refMeshVertices = inflator.vertices();
    refMeshElements = inflator.elements();



    if (dofOut != "") inflator.writePatternDoFs(dofOut + ".final.dof", result);

    cout << endl;
}

/* template<size_t _N, size_t _FEMDegree> */
/* void squareOptimizer (const po::variables_map &args, */ 
/* 		       const Job<_N> *job, */ 
/* 		       const int symmetryMode, */ 
/* 		       ETensor<_N>    & tensorC, */
/* 		       vector<double> & finalC, */ 
/* 		       vector<double> & finalParams) */
/* { */
/* 	vector<MeshIO::IOVertex> refMeshVertices; */
/* 	vector<MeshIO::IOElement> refMeshElements; */
  
/* 	Eigen::Matrix2d def; */

/* 	def << 1.0, 0.0, */
/* 		   0.0, 1.0; */

/* 	vector<double> dummyVec; */
/*     string solver = "levenberg_marquardt"; */
/* 	execute<_N, _FEMDegree>(args, job, def, symmetryMode, solver, refMeshVertices, refMeshElements, dummyVec, finalParams, true); */
	
/* 	// compute the effective C of the deformed cell */
/* 	auto &mat = HMG<_N>::material; */
/*     if (args.count("material")) mat.setFromFile(args["material"].as<string>()); */
/*     typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh; */
/*     typedef LinearElasticity::Simulator<Mesh> Simulator; */
/*     typedef typename Simulator::VField VField; */
/*     Simulator sim(refMeshElements, refMeshVertices); */
/*     auto &mesh = sim.mesh(); */

/*     vector<VField> w_ij; */
/*     mat.setTensor(mat.getTensor().transform(def.inverse())); */
/*     solveCellProblems(w_ij, sim); */
/*     auto effectiveC = homogenizedElasticityTensor(w_ij, sim).transform(def); */

/*     tensorC = effectiveC; */

/* 	double E1, E2, nu, mu; */
/* 	effectiveC.getOrthotropic2D(E1, E2, nu, mu); */
/* 	auto diffC = effectiveC - job->targetMaterial.getTensor(); */
/* 	auto diffS = effectiveC.inverse() - job->targetMaterial.getTensor().inverse(); */
/* 	double fNormC = diffC.quadrupleContract(diffC) / 2.0; */
/* 	double fNormS = diffS.quadrupleContract(diffS) / 2.0; */
	
/* 	finalC.clear(); */
/* 	finalC = {E1, E2, nu, mu, fNormC, fNormS}; */

/* 	vector<ETensor<_N>> tensorC_vec; */
/* 	vector<vector<double>> finalC_vec; */
	
/* 	tensorC_vec.push_back(tensorC); */
/* 	finalC_vec.push_back(finalC); */

	
/*  	// write the square meshes */
/* 	std::string refMeshName = "square_ref.msh"; */

/*     auto savePath = current_path() / args["out_folder"].as<string>() / "meshes"; */


/*     if (!exists(savePath) ||  !is_directory(savePath)) */
/*     	create_directory(savePath); */

/* 	MeshIO::save(savePath.string() + "/" + refMeshName, refMeshVertices, refMeshElements); */

/* 	// write the tensors and errors */
/* 	tableToFile(finalC_vec, savePath/"..", "sqr_optimization_results.txt"); */    
/* 	tableToFile(tensorC_vec, savePath/"..", "sqr_final_C.txt"); */    

/* } */

template<size_t _N, size_t _FEMDegree>
void cellPass (const po::variables_map &args, 
		       const Job<_N> *job, 
		       const Eigen::Matrix2d def, 
		       const double angle,
		       const int passID,
		       const int quadID,
		       vector<double> & outParams,
		       ETensor<_N>    & tensorC,
		       vector<double> & finalC)
{
	vector<MeshIO::IOVertex> refMeshVertices;
	vector<MeshIO::IOElement> refMeshElements;
 	 
	execute<_N, _FEMDegree>(args, job, def, quadID, refMeshVertices, refMeshElements, outParams);
	
	// compute the effective C of the deformed cell
	auto &mat = HMG<_N>::material;
    if (args.count("material")) mat.setFromFile(args["material"].as<string>());
    typedef LinearElasticity::Mesh<_N, _FEMDegree, HMG> Mesh;
    typedef LinearElasticity::Simulator<Mesh> Simulator;
    typedef typename Simulator::VField VField;
    Simulator sim(refMeshElements, refMeshVertices);
    auto &mesh = sim.mesh();

    vector<VField> w_ij;
    mat.setTensor(mat.getTensor().transform(def.inverse()));
    solveCellProblems(w_ij, sim);
    auto effectiveC = homogenizedElasticityTensor(w_ij, sim).transform(def);

    tensorC = effectiveC;

	double E1, E2, nu, mu;
	effectiveC.getOrthotropic2D(E1, E2, nu, mu);
	auto diffC = effectiveC - job->targetMaterial.getTensor();
	auto diffS = effectiveC.inverse() - job->targetMaterial.getTensor().inverse();
	double fNormC = diffC.quadrupleContract(diffC) / 2.0;
	double fNormS = diffS.quadrupleContract(diffS) / 2.0;
	
	finalC.clear();
	finalC = {E1, E2, nu, mu, fNormC, fNormS};

	// generate the deformed mesh
	auto bbox = mesh.boundingBox();
    VectorND<_N> center = 0.5 * (bbox.minCorner + bbox.maxCorner);
    vector<MeshIO::IOVertex> defMeshVertices;


	/* Eigen::Matrix2d rot, xReflect; */
	/* rot(0, 0) =  std::cos(angle); */
	/* rot(0, 1) = -std::sin(angle); */
	/* rot(1, 0) =  std::sin(angle); */
	/* rot(1, 1) =  std::cos(angle); */

	/* xReflect(0, 0) = -1.0; */
	/* xReflect(0, 1) =  0.0; */
	/* xReflect(1, 0) =  0.0; */
	/* xReflect(1, 1) =  1.0; */

	/* Eigen::Matrix2d transformation = rot * xReflect * def; */

    for (size_t vi = 0; vi < mesh.numVertices(); ++vi) {
        VectorND<_N> p = mesh.vertex(vi).node()->p;
        //defMeshVertices.emplace_back((transformation * (p - center) + center).eval());
        defMeshVertices.emplace_back((def * p).eval());
    }		

	// write the output meshes
	std::string defMeshName = "quad_" + to_string(quadID) + "_pass_" + to_string(passID) + "_def.msh";
	std::string refMeshName = "quad_" + to_string(quadID) + "_pass_" + to_string(passID) + "_ref.msh";

    auto savePath = current_path() / args["output"].as<string>() / "final_meshes";


    if (!exists(savePath) ||  !is_directory(savePath))
    	create_directory(savePath);

	MeshIO::save(savePath.string() + "/" + refMeshName, refMeshVertices, refMeshElements);
	MeshIO::save(savePath.string() + "/" + defMeshName, defMeshVertices, refMeshElements);
}

template<size_t _N, size_t _FEMDegree>
void meshPass (const po::variables_map &args, 
		       const Job<_N> *job, 
		       const vector<Eigen::Matrix2d> defs, 
		       const vector<double> angles,
		       const int passID,
		       vector<ETensor<_N>>    & tensorTable,
		       vector<vector<double>> & stiffnessTable,  
		       vector<vector<double>> & parameterTable)
{

	// find the deformation that conform to the symmetryMode best
	vector<Eigen::Matrix2d> newDefs;
	newDefs.clear();
	for (auto def = defs.begin(); def != defs.end(); ++def)
	{
		// extraxt relevant part of *def
		newDefs.push_back(*def);
	}


	stiffnessTable.clear();
	parameterTable.clear();
	tensorTable.clear();

	vector<double> parameterRow, stiffnessRow;
	ETensor<_N> stiffnessTensor;


	/* // run square optimization ... */
	/* squareOptimizer<_N, _FEMDegree>(args, job, symmetryMode, stiffnessTensor, stiffnessRow, squarParams); */

	int quadID = 0;
	for (auto def = newDefs.begin(); def != newDefs.end(); ++def)
	{
		// TODO: 
		// 1) setup a simulation with *def, 
		// 2) run it 
		// 3) store the in targetC and targetParams
		// 4) create the output percell meshes (prepend the passID to the file names)
		// 5) create the output tiled mesh (prepend the passID to the file name)
		cellPass<_N, _FEMDegree>(args, job, *def, angles[quadID], passID, quadID, parameterRow, stiffnessTensor, stiffnessRow);

		stiffnessTable.push_back(stiffnessRow);
		parameterTable.push_back(parameterRow);
		tensorTable.push_back(stiffnessTensor);
		++quadID;
	}
}

template<size_t _N, size_t _FEMDegree>
void tileQuad (const po::variables_map &args, 
		       const Job<_N> *job, 
		       const bool avgThick,
		       const int symmetryMode,
		       const vector<vector<double>> & parameterTable,
		       PolyMesh & pmesh)
{
	// tile the quad mesh using the parameterTable ... 
    //auto wi = WireInflator2D::construct(args["pattern"].as<string>());
    auto wi = WireInflator2D::construct<WireMesh2DMorteza>(args["pattern"].as<string>(), symmetryMode);
    WireInflator2D::OutMeshType mesh;

    CellParameters p_params = wi->createParameters();
    /* auto ops = wi->getParameterOperations(); */
    /* for (size_t i = 0; i < p_params.numberOfParameters(); ++i) { */
    /*     if (ops[i].type == ParameterOperation::Radius) */
    /*         range.push_back({minRadius, maxRadius}); */
    /*     else range.push_back(wi->getParameterRange(i)); */
    /* } */ 
    cout << parameterTable[0].size() << endl;
	cout << p_params.numberOfParameters() << endl;
    std::vector<CellParameters> quadParams;
    size_t nParams = p_params.numberOfParameters();
    for(size_t i = 0; i < parameterTable.size(); ++i)
	{
        assert(parameterTable[i].size() == nParams);
        quadParams.emplace_back(wi->numberOfParameters());
        auto &p = quadParams.back();
        for (size_t j = 0; j < nParams; ++j)
            p.parameter(j) = parameterTable[i][j];
    }
	
    TessellationParameters t_params;

    //static const bool averageThicknessOnBoundary = true;

    wi->generateQuadsPattern(pmesh, quadParams, t_params,
                             mesh, avgThick);

    std::vector<MeshIO::IOVertex> outVertices;
    std::vector<MeshIO::IOElement> outElements;
    for (const auto &p : mesh.nodes)
        outVertices.push_back(MeshIO::IOVertex(p[0], p[1], 0));
    for (const auto &e : mesh.elements)
        outElements.push_back(MeshIO::IOElement(e[0], e[1], e[2]));


	std::string tiledMeshName;
	if (avgThick)
		tiledMeshName = "tiled_avg_thick_yes.msh";
	else
		tiledMeshName = "tiled_avg_thick_no.msh";

	
    auto meshPath = current_path() / args["output"].as<string>() / "meshes";


    if (!exists(meshPath) ||  !is_directory(meshPath))
    	create_directory(meshPath);

	MeshIO::save(meshPath.string() + "/" + tiledMeshName, outVertices, outElements);

}



// helper functions to generate the output files/meshes
void tableToFile(const vector<ETensor<2>> & table,
                 path savePath,
                 const std::string fileName)
{
    ofstream out;
    out.open(savePath.string() + "/" + fileName);

    for (auto it = table.begin(); it != table.end(); ++it)
    {   
        out << *it << endl;
    }   
    out.close();
}

// helper functions to generate the output files/meshes
void tableToFile(const vector<vector<double>> & table,
				 path savePath,
				 const std::string fileName)
{
	ofstream out;
	out.open(savePath.string() + "/" + fileName);

	for (auto it = table.begin(); it != table.end(); ++it)
	{
		vector<double> row = *it;
		for (size_t i = 0; i < row.size() - 1; ++i)
		{
			out << showpos << scientific;
			out << row[i] << "\t";
		}
		out << row.back() << "\n";
	}
	out.close();
}

void tableToFile(const vector<Eigen::Matrix2d> & table,
				 path savePath,
				 const std::string fileName)
{

	ofstream out;
	out.open(savePath.string() + "/" + fileName);

	for (auto it = table.begin(); it != table.end(); ++it)
	{
		Eigen::Matrix2d mat = *it;
		out << showpos << scientific <<  mat(0, 0) << "\t"  
			<< showpos << scientific <<  mat(0, 1) << "\t" 
			<< showpos << scientific <<  mat(1, 0) << "\t" 
			<< showpos << scientific <<  mat(1, 1) << "\n";
	}
	out.close();
}


void writeTilerParameters(const vector<vector<double>> & stiffnessTable, 
		         		  const vector<vector<double>> & parameterTable,
		         		  path savePath,
		         		  const std::string fileName)
{
	assert(stiffnessTable.size() == parameterTable.size());

	vector<vector<double>> table;
	table.clear();

	for (size_t i = 0; i < stiffnessTable.size(); ++i)
	{
		vector<double> row;
		row.clear();
		for (int j = 0; j < 5; ++j)
			(j == 4) ? row.push_back(j + 5) : row.push_back(j);
		for (size_t j = 0; j < parameterTable[i].size(); ++j)
			row.push_back(parameterTable[i][j]);
		table.push_back(row);
	}
	tableToFile(table, savePath, fileName);
}

void checkPmesh(PolyMesh & pmesh)
{
	typedef WireMeshEmbedding<EMesh, PolyMesh>				WireEmbedding;
	typedef typename WireEmbedding::QuadParametrization		QuadParametrization;
	int faceNumber = 0;
	for (auto fc = pmesh.face.begin(); fc != pmesh.face.end(); ++fc)
	{
		QuadParametrization qpar = WireEmbedding::getQuadParametrizationHandle(pmesh)[fc];
		cout << "index0 for face " << faceNumber << " is " << int(qpar.index0) << endl;
 
		++faceNumber;
	}
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }   
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

vector<Eigen::Matrix2d>  truncate(const int decimalPlace, vector<Eigen::Matrix2d> mats)
{
	vector<Eigen::Matrix2d> outMats;
	outMats.clear();
	long int factor = pow(10, decimalPlace);
	for (size_t i = 0; i < mats.size(); ++i)
	{
		Eigen::Matrix2d mat = mats[i];
		Eigen::Matrix2d outMat;
		outMat(0, 0) = std::trunc(mat(0, 0) * factor) / factor;
		outMat(0, 1) = std::trunc(mat(0, 1) * factor) / factor;
		outMat(1, 0) = std::trunc(mat(1, 0) * factor) / factor;
		outMat(1, 1) = std::trunc(mat(1, 1) * factor) / factor;

		outMats.push_back(outMat);
	}
	return outMats;
}

/* vector<vector<double>> rowMajor(const vector<Eigen::Matrix2d> inDefs) */
/* { */
/* 	vector<vector<double>> outDefs; // in row major */
/* 	for (auto currentDef : inDefs) */
/* 	{ */
/* 		Eigen::Matrix<double, 2, 2, RowMajor> currentDefRowMajor = *currentDef; */
/* 		vector<double> currentRow; */
/* 		currentRow.clear(); */
/* 		for (size_t i = 0; i < currentDefRowMajor.size(); ++i) */
/* 			currentRow.push_back(currentDefRowMajor.data()[i]); */
/* 		outDefs.push_back(currentRow); */
/* 	} */
/* 	return outDefs; */
/* } */

vector<double> processDeformation(const Eigen::Matrix2d F)
{
	vector<double> decomposedF; // {lambda1, lambda2, theta, psi, such that R[psi]*Q[theta]*Diag(lambda1,lambda2)*Q^T[theta]=F

	/* cout << "F is: " << endl << F << endl; */
	/* // find U and R first: U^2=F^T*F  and  R=F*U^-1 */
	/* Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es1(F.transpose() * F); */
	
	/* Eigen::Matrix2d U = es1.operatorSqrt(); */
	/* cout << "U is: " << endl << U << endl; */
	/* Eigen::Matrix2d R = F * U.inverse(); */
	/* Eigen::Rotation2Dd R_rot; */
	/* R_rot.fromRotationMatrix(R); */
	/* cout << "angle of R is: " << R_rot.angle() * 180 / std::acos(-1.0) << endl; */
	/* cout << "R is: " << endl << R << endl; */
	/* cout << "R_rot is: " << endl << R_rot.matrix() << endl; */
	/* cout << "U is: " << endl << U << endl; */
	/* double psi = std::atan2(R(1, 0), R(0, 0)) * 180 / std::acos(-1.0); */

	/* // find Q and D: */ 
	/* Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> es2(U); */

	/* Eigen::Vector2d eVals = es2.eigenvalues(); */
	/* Eigen::Matrix2d     Q = es2.eigenvectors(); */
	/* Eigen::Rotation2D<double> Q_rot; */
	/* Q_rot.fromRotationMatrix(Q); */
	/* cout << "angle of Q is: " << Q_rot.angle() * 180 / std::acos(-1.0) << endl; */
	/* cout << "check is " << endl << Q.transpose() * U * Q << endl; */
	/* cout << "Q is: " << endl << Q << endl; */
	/* cout << "Q_rot is: " << endl << Q_rot.matrix() << endl; */
	/* cout << "|Q_rot - Q| is: " << (Q - Q_rot.matrix()).norm() << endl; */

	/* double theta = std::atan2(Q(1, 0), Q(0, 0)) * 180 / std::acos(-1.0); */

	/* decomposedF = {eVals(0), eVals(1), theta, psi}; */
	
	
	return decomposedF;
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    //cout << setprecision(16);
	ofstream ofs;
	ofs.open("test_output.txt");
	ofs << "grad\\_JS" << "\t" << "grad\\_RT" << "\t" << "JS" << "\t" << "RT" << std::endl;
	ofs.close();

	// read in the cmd aurguments
    po::variables_map args = parseCmdLine(argc, argv);
	
	// make the output directory
    auto savePath = current_path() / args["output"].as<string>();
    if (!exists(savePath) ||  !is_directory(savePath))
    	create_directory(savePath);
   
   	// read the job file
    auto job = parseJobFile(args["job"].as<string>());

    int symmetryMode = args["sym"].as<int>();



	vector<Eigen::Matrix2d> defs;
	defs.clear();
	vector<double> angles;
	Eigen::Matrix2d def;
	std::vector<double> lambdas = {1.0, 1.3, 1.5};
	
	// different stretches
	for (auto lam = lambdas.begin(); lam != lambdas.end(); ++ lam){
		def << *lam , 0.0,
			    0.0 , 1.0;
		defs.push_back(def);
		def << 1.0 , 0.0,
			   0.0 , *lam;
		defs.push_back(def);
		angles.push_back(0.0);
		angles.push_back(0.0);
	}

	tableToFile(defs, savePath, "defs_original.txt");

	tableToFile(truncate(3, defs), savePath, "defs_truncated.txt");
	vector<vector<double>> optimizationTable, parameterTable;
	vector<ETensor<2>> tensorTable;
	auto job2D = dynamic_cast<Job<2> *>(job);
	meshPass<2, 2>(args, job2D, 
			       truncate(3, defs), angles, 
			       1, 
			       tensorTable, optimizationTable, parameterTable);

	vector<vector<double>> optimizationTable1, optimizationTable2;
	vector<ETensor<2>> tensorTable1, tensorTable2;

	optimizationTable1.clear();
	optimizationTable2.clear();
	tensorTable1.clear();
	tensorTable2.clear();

	for (size_t i = 0; i < lambdas.size(); ++i)
	{
		optimizationTable1.push_back(optimizationTable[2*i]);
		optimizationTable2.push_back(optimizationTable[2*i+1]);

		tensorTable1.push_back(tensorTable[2*i]);
		tensorTable2.push_back(tensorTable[2*i+1]);
	}

	tableToFile(optimizationTable1, savePath, "gpr_01_optimization_results.txt");    
	tableToFile(optimizationTable2, savePath, "gpr_02_optimization_results.txt");    
	tableToFile(tensorTable1, savePath, "gpr_01_final_C.txt");    
	tableToFile(tensorTable2, savePath, "gpr_02_final_C.txt");    
    return 0;
}
