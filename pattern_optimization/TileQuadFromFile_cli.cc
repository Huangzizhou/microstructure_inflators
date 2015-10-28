////////////////////////////////////////////////////////////////////////////////
//// TileQuadFromFile.cc
//////////////////////////////////////////////////////////////////////////////////
///*! @file
////      This reads in a (2D) quad mesh file and a txt file containing the
////      optimized pattern params and generates the tiled mesh
//*/ 
////  Author:  Morteza H Siboni , hakimi1364@gmail.com
////  Note:    This is written based on Julian's PatternOptimization_cli
////  Company:  New York University
////  Created:  OCT 2015 
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
    cout << "Usage: TileQuadFromFile_cli [options] quadMesh.obj" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}


po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("quad", po::value<string>(), "quadMesh.obj file");

    po::positional_options_description p;
    p.add("quad", 1);

    po::options_description visible_opts;
   	visible_opts.add_options()("help",			"Produce this help message")
        ("pattern,p",			po::value<string>(), "Pattern wire mesh (.obj|wire)")
        ("material,m",			po::value<string>(), "base material")
        ("inParam",				po::value<string>(), "input optimized parameters file")
        ("inStiff",				po::value<string>()->default_value(""), "input optimized stiffnesses file")
        ("inCost",	 			po::value<string>()->default_value(""), "input costs file")
        ("output,o",     		po::value<string>(), "final tiled .msh file")
        ("max_volume,v", 		po::value<double>(), "maximum element volume parameter for wire inflator")
        ("avg_thickness,t", 	po::value<bool>()->default_value(true), "true to average the thicknesses when tiling")
        ("sym",          		po::value<int>()->default_value(3), "symmetry mode")
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
		cout << "Error: must specify a quadMesh.obj" << endl;
		fail = true;
	}

	if (vm.count("pattern") == 0) {
		cout << "Error: must specify the edge mesh!" << endl;
		fail = true;
	}

	if (vm.count("output") == 0) {
		cout << "Error: must specify output.msh file" << endl;
		fail = true;
	}

	if (fail || vm.count("help"))
		usage(fail, visible_opts);

	return vm;
}

template<size_t cols>
void readTable(std::string fileName, vector<vector<Real>> & dataTable)
{
	cout << endl << "# of cols is" << cols << endl;
	ifstream in;
	in.open(fileName);
	string line;
	
	while (getline(in, line)){
		
		vector<string> numbers;
		boost::trim(line);
		boost::split(numbers, line, boost::is_any_of("\t "),
					 boost::token_compress_on);



		if (numbers.size() != cols)
			throw runtime_error("invalid number of columns");
		
		vector<Real> row;
		for (size_t i = 0; i < numbers.size(); ++i)
			row.push_back(stod(numbers[i]));
		
		dataTable.push_back(row);
	}

	cout << "read " << dataTable.size() << " lines from " << fileName << "." << endl << endl;

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

template<size_t _N>
void tileQuad (const po::variables_map &args)
{
	
	// read in the quadMesh
    PolyMesh pmesh;
    readQuadMesh(args["quad"].as<string>(), pmesh);
	// tile the quad mesh using the parameterTable ... 
    auto wi = WireInflator2D::construct(args["pattern"].as<string>(), args["sym"].as<int>());
    WireInflator2D::OutMeshType mesh;

    CellParameters p_params = wi->createParameters();
 
 	std::vector<CellParameters> quadParams;
    size_t nParams = p_params.numberOfParameters();
	
	vector<vector<Real>> parameterTable;
	switch (nParams){
		case 9:
			readTable<9>(args["inParam"].as<string>(), parameterTable);
			break;
		case 13:
			break;
	}


    for(size_t i = 0; i < parameterTable.size(); ++i)
	{
        assert(parameterTable[i].size() == nParams);
        quadParams.emplace_back(wi->numberOfParameters());
        auto &p = quadParams.back();
        for (size_t j = 0; j < nParams; ++j)
            p.parameter(j) = parameterTable[i][j];
    }
	
    TessellationParameters t_params;
    t_params.max_area = args["max_volume"].as<double>();
    //static const bool averageThicknessOnBoundary = true;

    wi->generateQuadsPattern(pmesh, quadParams, t_params,
                             mesh, args["avg_thickness"].as<bool>());

    std::vector<MeshIO::IOVertex> outVertices;
    std::vector<MeshIO::IOElement> outElements;
    for (const auto &p : mesh.nodes)
        outVertices.push_back(MeshIO::IOVertex(p[0], p[1], 0));
    for (const auto &e : mesh.elements)
        outElements.push_back(MeshIO::IOElement(e[0], e[1], e[2]));


	std::string tiledMeshName = args["output"].as<string>();
	MeshIO::save(tiledMeshName, outVertices, outElements);

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
	// read in the cmd aurguments
    po::variables_map args = parseCmdLine(argc, argv);
    tileQuad<2>(args);

    return 0;
}
