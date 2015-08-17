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
        ("out_folder,f",				po::value<string>(), "folder name to store the results")
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
	
	if (vm.count("out_folder") == 0) {
		cout << "Error: must specify a name for the output folder!" << endl;
		fail = true;
	}

	if (fail || vm.count("help"))
		usage(fail, visible_opts);

	return vm;
}


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
			 const int symMode, 
			 vector<MeshIO::IOVertex> & refMeshVertices,
             vector<MeshIO::IOElement> &  refMeshElements,
			 vector<double> & finalParams) 
{
    shared_ptr<ConstrainedInflator<_N>> inflator_ptr;
    if (_N == 2) {
        /* inflator_ptr = make_shared<ConstrainedInflator<_N>>( */
        /*         job->parameterConstraints, */
        /*         args["pattern"].as<string>(), */
        /*         args["sym"].as<int>()); // the last parameter is a symmetryMode (see the corresponding constructor in Inflator.hh for more details) */
        
        // this uses the original Luigi's symmetry stuff
        inflator_ptr = make_shared<ConstrainedInflator<_N>>(
                job->parameterConstraints,
                args["pattern"].as<string>()); // MHS JUL14, 2015: the last parameter is a symmetryMode (see the corresponding constructor in Inflator.hh for more details)
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
    mat.setTensor(mat.getTensor().transform(deformation.inverse()));

    string dofOut = args["dofOut"].as<string>();
    if (dofOut != "")
        inflator.setDoFOutputPrefix(dofOut);

   	SField params(job->initialParams);
   	
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

    //std::cout << "Final p:";
    std::vector<Real> result(params.domainSize());
    finalParams.clear();
    for (size_t i = 0; i < result.size(); ++i) {
        result[i] = params[i];
        //cout << "\t" << params[i];
        finalParams.push_back(params[i]);
    }

    refMeshVertices = inflator.vertices();
    refMeshElements = inflator.elements();



    if (dofOut != "") inflator.writePatternDoFs(dofOut + ".final.dof", result);

    cout << endl;
}


template<size_t _N, size_t _FEMDegree>
void cellPass (const po::variables_map &args, 
		       const Job<_N> *job, 
		       const Eigen::Matrix2d def, 
		       const double angle,
		       const int symmetryMode, 
		       const int passID,
		       const int quadID,
		       vector<double> & finalC, 
		       vector<double> & finalParams)
{
	vector<MeshIO::IOVertex> refMeshVertices;
	vector<MeshIO::IOElement> refMeshElements;
  
	execute<_N, _FEMDegree>(args, job, def, symmetryMode, refMeshVertices, refMeshElements, finalParams);
	
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

    auto savePath = current_path() / args["out_folder"].as<string>() / "meshes";


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
		       const int symmetryMode,
		       const int passID,
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

	vector<double> stiffnessRow, parameterRow;
	
	int quadID = 0;
	for (auto def = newDefs.begin(); def != newDefs.end(); ++def)
	{
		// TODO: 
		// 1) setup a simulation with *def, 
		// 2) run it 
		// 3) store the in targetC and targetParams
		// 4) create the output percell meshes (prepend the passID to the file names)
		// 5) create the output tiled mesh (prepend the passID to the file name)
		cellPass<_N, _FEMDegree>(args, job, *def, angles[quadID], symmetryMode, passID, quadID, stiffnessRow, parameterRow);

		stiffnessTable.push_back(stiffnessRow);
		parameterTable.push_back(parameterRow);
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
    auto wi = WireInflator2D::construct(args["pattern"].as<string>());
    WireInflator2D::OutMeshType mesh;

    CellParameters p_params = wi->createParameters();
    /* auto ops = wi->getParameterOperations(); */
    /* for (size_t i = 0; i < p_params.numberOfParameters(); ++i) { */
    /*     if (ops[i].type == ParameterOperation::Radius) */
    /*         range.push_back({minRadius, maxRadius}); */
    /*     else range.push_back(wi->getParameterRange(i)); */
    /* } */ 

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

	
    auto meshPath = current_path() / args["out_folder"].as<string>() / "meshes";


    if (!exists(meshPath) ||  !is_directory(meshPath))
    	create_directory(meshPath);

	MeshIO::save(meshPath.string() + "/" + tiledMeshName, outVertices, outElements);

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
			out << row[i] << "\t";
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
		out << mat(0, 0) << "\t" << mat(0, 1) << "\t" <<mat(1, 0) << "\t" <<mat(1, 1) << "\n";
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

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    //cout << setprecision(16);


	// read in the cmd aurguments
    po::variables_map args = parseCmdLine(argc, argv);


    auto job = parseJobFile(args["job"].as<string>());

    auto savePath = current_path() / args["out_folder"].as<string>();


    if (!exists(savePath) ||  !is_directory(savePath))
    	create_directory(savePath);


	std::string outPathStr = args["out_folder"].as<string>();
    

	// read in the quadMesh
    PolyMesh pmesh;
    readQuadMesh(args["quad"].as<string>(), pmesh);

	// compute the stretches/jacobians for each quad in the quadMesh
	vector<Eigen::Matrix2d> defs;
	vector<double>			angles;
	getDeformations(pmesh, defs, angles, "jacobian");


	for (size_t i = 0; i < angles.size(); ++i)
		cout << "angle of quad " << i+1 << " is " << angles[i] * 180.0 / 3.141595  << endl;

	tableToFile(defs, savePath, "defs.txt");

	vector<vector<double>> stiffnessTable, parameterTable;

	auto job2D = dynamic_cast<Job<2> *>(job);

	meshPass<2, 2>(args, job2D, defs, angles, 1, 0, stiffnessTable, parameterTable);


	writeTilerParameters(stiffnessTable, parameterTable, savePath, "tile_params.txt");

	tileQuad<2, 2> (args, job2D, true,  1, parameterTable, pmesh);
	tileQuad<2, 2> (args, job2D, false, 1, parameterTable, pmesh);

    return 0;
}
