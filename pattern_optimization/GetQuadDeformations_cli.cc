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
/* #include <Materials.hh> */
/* #include <PeriodicHomogenization.hh> */
#include "GlobalBenchmark.hh"
#include "PolyMeshType.h"
#include "PolyMeshUtils.h"
#include "WireMeshEmbedding.h"

#include "EdgeMeshType.h"
/* #include "Inflator.hh" */

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


/* #include "PatternOptimization.hh" */
/* #include "PatternOptimizationJob.hh" */

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace fs;
using namespace std;
/* using namespace PatternOptimization; */
/* using namespace PeriodicHomogenization; */

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: GetQuadDeformations_cli [options] quadMesh.obj" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}


po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("quad", po::value<string>(), "input quad mesh file");

    po::positional_options_description p;
    p.add("quad", 1);

    po::options_description visible_opts;
   	visible_opts.add_options()("help",			"Produce this help message")
        ("output,o",     po::value<string>()->default_value(""),                 "output .txt of deformations F11 F12 F21 F22 per line")
        ("mode,m",       po::value<string>()->default_value("jacobian"),         "mode = {jacobian, stretch} for outputing the jacobian or right stretch tensor")
        ("truncate,t",   po::value<size_t>()->default_value(6),                  "truncates components of the deformation to the t decimal place")
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
	if (vm.count("quad") == 0) {
		cout << "Error: must specify input quadMesh.obj file!" << endl;
		fail = true;
	}

	if (vm.count("output") == 0) {
		cout << "Error: must specify the output.txt file!" << endl;
		fail = true;
	}
		
	if (fail || vm.count("help"))
		usage(fail, visible_opts);

	return vm;
}

template<size_t _N> 
using ETensor = ElasticityTensor<Real, _N>; 

typedef ScalarField<Real> SField;

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
		             const string type = "jacobian")
{
	typedef WireMeshEmbedding<EMesh, PolyMesh>				WireEmbedding;

	std::vector<double> angles;

	if (type == "jacobian")
		WireEmbedding::getJacobians(pmesh, deformations, angles);
	else if (type == "stretch")
		WireEmbedding::getStretches(pmesh, deformations, angles);
	else
		;
}

vector<Eigen::Matrix2d> truncate(const int decimalPlace, vector<Eigen::Matrix2d> mats) 
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

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
	// read in the cmd aurguments
    po::variables_map args = parseCmdLine(argc, argv);

    auto savePath = current_path();

	std::string outPathStr = args["output"].as<string>();
    

	// read in the quadMesh
    PolyMesh pmesh;
    readQuadMesh(args["quad"].as<string>(), pmesh);

	// compute the stretches/jacobians for each quad in the quadMesh
	vector<Eigen::Matrix2d> defs;
	getDeformations(pmesh, defs, args["mode"].as<string>());

	// truncate the stretches/jacobians
	defs = truncate(args["truncate"].as<size_t>(), defs);	
		

	tableToFile(defs, savePath, args["output"].as<string>());

    return 0;
}
