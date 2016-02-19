////////////////////////////////////////////////////////////////////////////////
////  GetQuadDeformations_cli.cc
//////////////////////////////////////////////////////////////////////////////////
///*! @file
////      This reads in a (2D) quad mesh, finds the equivalent parallelograms, and
////	  returns the deformation (or stretch) of each parallelogram 
//*/ 
////  Author  : Morteza H Siboni (mhs), m.hakimi.siboni.@gmail.com
////  Company : New York University
////  Created : 08/10/2014
//////////////////////////////////////////////////////////////////////////////////
#include "PolyMeshType.h"
#include "PolyMeshUtils.h"
#include "WireMeshEmbedding.h"
#include "EdgeMeshType.h"
#include "table.hh"

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <algorithm>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace fs;
using namespace std;

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
        ("allDefs,a", 	      po::value<string>(),								"a .txt: including all deformations F11 F12 F21 F22")
        ("uniqueDefs,u",      po::value<string>(),								"a .txt: including unique copy of all deformations")
        ("conversionMap,c",   po::value<string>(), 								"a .txt: including the conversion map from all-defs to unique-defs") 
        ("mode,m",       po::value<string>()->default_value("jacobian"),	"mode = {jacobian, stretch}")
        ("truncate,t",   po::value<size_t>()->default_value(4),				"truncates components of the deformation to the t decimal place")
        ("matlabData,d", po::value<string>()->default_value(""), 			"output data file name for matlab visualization")
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

	if (vm.count("allDefs") == 0) {
		cout << "Error: must specify the .txt file to store all deformations!" << endl;
		fail = true;
	}

	if (vm.count("uniqueDefs") == 0) {
		cout << "Error: must specify the .txt file to store unique deformations!" << endl;
		fail = true;
	}	
	
	if (vm.count("conversionMap") == 0) {
		cout << "Error: must specify the .txt file to store the conversion map!" << endl;
		fail = true;
	}	
	if (fail || vm.count("help"))
		usage(fail, visible_opts);

	return vm;
}

typedef WireMeshEmbedding<EMesh, PolyMesh>				WireEmbedding;
typedef typename WireEmbedding::QuadParametrization 	QuadParametrization;


void dumpQuadCoord2Stream(PolyMesh & pmesh, bool parametrized=true, ostream &os = cout)
{
	int firstIdx;

	for(auto fc = pmesh.face.begin(); fc != pmesh.face.end(); ++fc)
	{
		QuadParametrization qp = WireEmbedding::getQuadParametrizationHandle(pmesh)[fc];

		parametrized ? firstIdx = qp.index0 : firstIdx = 0;

		for (size_t i = 0; i < 4; ++i)
		{
			int updatedIdx = (i + firstIdx) % 4;
			os << fc->cP(updatedIdx)[0] << "\t" << fc->cP(updatedIdx)[1] << endl;
		}
	}
}

// read the quad mesh:
void readQuadMesh(const string & quadMeshPath, PolyMesh & pmesh)
{
	typedef PolyMeshUtils<PolyMesh>		PMU;
	bool ok = false;
	ok = PMU::importFromOBJ(quadMeshPath, pmesh);
	if (ok)
		WireEmbedding::preprocessQuadMesh(pmesh);

	/* WireEmbedding::dumpParametrizationSequence(pmesh); // this is meant for debugging purposes */ 
	WireEmbedding::createLocalParametrization(pmesh); // overwrite Luigi's coherent parametrization 
	/* WireEmbedding::dumpParametrizationSequence(pmesh); // this is meant for debugging purposes */ 
}

// get the deformation (F or U) for each quad in the quadMesh 
void getDeformations(PolyMesh & pmesh, 
		             vector<Eigen::Matrix2d> & deformations, 
		             const string type = "jacobian")
{
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


/* void tableToFile(const vector<Eigen::Matrix2d> & table, */
/* 				 path savePath, */
/* 				 const std::string fileName) */
/* { */

/* 	ofstream out; */
/* 	out.open(savePath.string() + "/" + fileName); */

/* 	for (auto it = table.begin(); it != table.end(); ++it) */
/* 	{ */
/* 		Eigen::Matrix2d mat = *it; */
/* 		out << showpos << scientific <<  mat(0, 0) << "\t" */
/* 			<< showpos << scientific <<  mat(0, 1) << "\t" */
/* 			<< showpos << scientific <<  mat(1, 0) << "\t" */
/* 			<< showpos << scientific <<  mat(1, 1) << "\n"; */
/* 	} */
/* 	out.close(); */
/* } */

/* void tableToFile(const vector<vector<double>> & table, */
/* 				 path savePath, */
/* 				 const std::string fileName) */
/* { */

/* 	ofstream out; */
/* 	out.open(savePath.string() + "/" + fileName); */

/* 	for (auto it = table.begin(); it != table.end(); ++it) */
/* 	{ */
/* 		vector<double> row = *it; */
/* 		for (size_t i = 0; i < row.size() - 1; ++i) */
/* 			out << showpos << scientific <<  row[i] << "\t"; */
/* 		out << showpos << scientific <<  row[row.size() - 1] << "\n"; */
/* 	} */
/* 	out.close(); */
/* } */


void dumpTable(const vector<vector<double>> & table, ostream & os)
{

	for (auto it = table.begin(); it != table.end(); ++it)
	{
		vector<double> row = *it;
		for (size_t i = 0; i < row.size() - 1; ++i)
			os << showpos << scientific <<  row[i] << "\t";
		os << showpos << scientific <<  row[row.size() - 1] << "\n";
	}
}

void eigenTable2vecTable(const vector<Eigen::Matrix2d> eigenTable, vector<vector<double>> & vecTable)
{
	for (auto it = eigenTable.begin(); it != eigenTable.end(); ++it)
	{
		Eigen::Matrix2d mat = *it;
		vector<double> row;
		row.push_back(mat(0, 0));
		row.push_back(mat(0, 1));
		row.push_back(mat(1, 0));
		row.push_back(mat(1, 1));

		vecTable.push_back(row);
	}
}

void dumpMap(const map<int, int> & inMap, ostream & os)
{
	for (auto & p : inMap)
		os << p.first << "\t" << p.second << endl;
}

vector<vector<double>> checkMap(map<int, int> inMap, vector<vector<double>> longTable, vector<vector<double>> shortTable)
{
	vector<vector<double>> outTable;
	outTable = longTable;
	for (auto & p : inMap)
	{
		int idInLong  = p.first;
		int idInShort = p.second;

		for (size_t i = 0; i < outTable[idInLong].size(); ++i)
			outTable[idInLong][i] = outTable[idInLong][i] - shortTable[idInShort][i]; 
	}
	return outTable;
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

	// read in the quadMesh
    PolyMesh pmesh;
    readQuadMesh(args["quad"].as<string>(), pmesh);

    // output quad coordinates for visualization purposes 
    if (args["matlabData"].as<string>() != "" && args["matlabData"].as<string>() != "stdout"){
    	ofstream out_original;
    	ofstream out_parametrized;
    	
    	out_original.open(args["matlabData"].as<string>() + "_original.txt");
    	out_parametrized.open(args["matlabData"].as<string>() + "_parametrized.txt");
    
    	dumpQuadCoord2Stream(pmesh , false , out_original);
    	dumpQuadCoord2Stream(pmesh , true  , out_parametrized);

    	out_original.close();
    	out_parametrized.close();
	}
	else if (args["matlabData"].as<string>() == "stdout"){
		cout << "coordinates (x,y) for each quad in the original quad mesh:   " << endl;
		cout << "lines (i-1)*4+j, j=[1..4] corresponds to the ith quad's (x,y)" << endl;
		cout << "-------------------------------------------------------------" << endl;
		dumpQuadCoord2Stream(pmesh, false);
		cout << endl;

		cout << "coordinates (x,y) for each quad in the parametrized quad mesh:" << endl;
		cout << "lines (i-1)*4+j, j=[1..4] corresponds to the ith quad's (x,y) " << endl;
		cout << "--------------------------------------------------------------" << endl;
		dumpQuadCoord2Stream(pmesh, true);
		cout << endl;
	}
	else
		;


	// compute the stretches/jacobians for each quad in the quadMesh
	vector<Eigen::Matrix2d> defs;
	getDeformations(pmesh, defs, args["mode"].as<string>());

	// truncate the stretches/jacobians
	defs = truncate(args["truncate"].as<size_t>(), defs);	

	table<double, 2> defTable(defs);



	// getting/writing all the quad deformations to file
	vector<vector<double>> allDefs    = defTable.getTable();
	ofstream out;
	out.open(args["allDefs"].as<string>());
	dumpTable(allDefs, out);
	out.close();

	// getting/writing unique quad deformations to file
	vector<vector<double>> uniqueDefs = defTable.getUniqueTable();
	out.open(args["uniqueDefs"].as<string>());
	dumpTable(uniqueDefs, out);
	out.close();

	// getting/writing the conversion map to file
	map<int, int> conversionMap       = defTable.getMap();
	out.open(args["conversionMap"].as<string>());
	dumpMap(conversionMap, out);
	out.close();

    return 0;
}
