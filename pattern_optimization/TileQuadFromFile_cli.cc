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
#include <MSHFieldWriter.hh>
#include <MSHFieldParser.hh>
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
#include <boost/algorithm/string.hpp> // required for parsing jacobian

// the following are not required anymore
#include "PatternOptimization.hh" 
#include "PatternOptimizationJob.hh"

#include <vcg/complex/complex.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace fs;
using namespace std;

// the following are not required anymore
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
        ("scale,s",			 	po::value<double>()->default_value(1.0), "scales the input mesh by a factor of 's' [used when there is a significant change of scale in the quad mesh]")
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

void readTable(std::string fileName, vector<vector<Real>> & dataTable)
{
	ifstream in;
	in.open(fileName);
	string line;
	
	while (getline(in, line)){
		
		vector<string> numbers;
		boost::trim(line);
		boost::split(numbers, line, boost::is_any_of("\t "),
					 boost::token_compress_on);

		cout << numbers.size() << endl;	
		vector<Real> row;
		for (size_t i = 0; i < numbers.size(); ++i)
			row.push_back(stod(numbers[i]));
		
		dataTable.push_back(row);
	}
	cout << "----------" << endl;
	cout << "read " << dataTable.size() << " lines from " << fileName << "." << endl;
	cout << "----------" << endl;
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

	/* int faceCounter = 0; */
	/* for(auto fc = pmesh.face.begin(); fc != pmesh.face.end(); ++fc) */
	/* { */
	/* 	cout << "points in face " << faceCounter << " are: " << "("  << fc->cP(0)[0] << "," << fc->cP(0)[1] << ") " << */
	/* 		                                                    "("  << fc->cP(1)[0] << "," << fc->cP(1)[1] << ") " << */
	/* 		                                                    "("  << fc->cP(2)[0] << "," << fc->cP(2)[1] << ") " << */
	/* 		                                                    "("  << fc->cP(3)[0] << "," << fc->cP(3)[1] << ")"  << endl; */ 
	/* 	++faceCounter; */
	/* } */
}

// scale the quad mesh:
void scaleQuadMesh(double scaleFactor, PolyMesh & pmesh)
{
	for (auto vert = pmesh.vert.begin(); vert != pmesh.vert.end(); ++vert)
		vert->P() *= scaleFactor;
}

typedef typename PolyMesh::CoordType        PCoordType;
typedef typename PolyMesh::FaceType        	PFaceType;

Real quadArea(PolyMesh & pmesh, PFaceType & f)
{
	
	/* typedef WireMeshEmbedding<EMesh, PolyMesh>				WireEmbedding; */
	/* typedef typename WireEmbedding::QuadParametrization		QuadParametrization; */

	
	/* QuadParametrization qpar = getQuadParametrizationHandle(pmesh)[f]; */

	char index0 = 0;
	char index1 = (index0 + 1) % 4; 
	char index2 = (index0 + 2) % 4; 
	char index3 = (index0 + 3) % 4; 

	PCoordType p0 = f.cP(index0);
	PCoordType p1 = f.cP(index1);
	PCoordType p2 = f.cP(index2);
	PCoordType p3 = f.cP(index3);

	PCoordType a = p1 - p0;
	PCoordType b = p2 - p1;
	PCoordType c = p3 - p2;
	PCoordType d = p0 - p3;

	PCoordType m = p3 - p1;
	PCoordType n = p2 - p0;

	Real a2 = a[0] * a[0] + a[1] * a[1];
	Real b2 = b[0] * b[0] + b[1] * b[1];
	Real c2 = c[0] * c[0] + c[1] * c[1];
	Real d2 = d[0] * d[0] + d[1] * d[1];
          
	Real m2 = m[0] * m[0] + m[1] * m[1];
	Real n2 = n[0] * n[0] + n[1] * n[1];

	Real x = b2 + d2 - a2 - c2;

	return std::sqrt(4.0 * m2 * n2 - x * x ) / 4.0;
}



void findQuadIdsForElements(const std::vector<MeshIO::IOVertex>  &verts,
		                    const std::vector<MeshIO::IOElement> &elems,
		                    PolyMesh &pmesh, 
		                    std::vector<size_t> & ids,
		                    std::vector<Real>   & vols)
{

	typedef ClipperLib::IntPoint 				IntPoint;
	

	const int multiplier = (1 << 24);

	std::vector<IntPoint> elementCenters;

	for (size_t i = 0; i < elems.size(); ++i)
	{
		int idx0 = elems[i][0];
		int idx1 = elems[i][1];
		int idx2 = elems[i][2];
		Point3D p0 = verts[idx0];
		Point3D p1 = verts[idx1];
		Point3D p2 = verts[idx2];

		double centerX = (p0[0] + p1[0] + p2[0]) / 3.0;
		double centerY = (p0[1] + p1[1] + p2[1]) / 3.0;

		elementCenters.push_back(ClipperLib::IntPoint(ClipperLib::cInt(multiplier * centerX), ClipperLib::cInt(multiplier * centerY)));
	}

	std::vector<ClipperLib::Path> facePaths;
	std::vector<Real> faceVolumes;
	for (auto f = pmesh.face.begin(); f != pmesh.face.end(); ++f)
	{
		// get the clipping quad
		ClipperLib::Path facePath;
		for (char i=0; i<f->VN(); i++)
		{
			IntPoint pt = ClipperLib::IntPoint(ClipperLib::cInt(multiplier * f->cP(i)[0]), ClipperLib::cInt(multiplier * f->cP(i)[1]));
			facePath.push_back(pt);

		}

		facePaths.push_back(facePath);
		faceVolumes.push_back(quadArea(pmesh, *f));
	}

	for (size_t i = 0; i < elementCenters.size(); ++i)
		for (size_t j = 0; j < facePaths.size(); ++j)
		{
			int isInside = ClipperLib::PointInPolygon(elementCenters[i], facePaths[j]);
			int flag = 0;
			if (isInside != 0)
			{
				flag = 1;
				ids.push_back(j);
				vols.push_back(faceVolumes[j]);
			}
			if (flag == 1)
				break;
		}

	Real totalVolume = 0;
	for (size_t i = 0; i < faceVolumes.size(); ++i)
		totalVolume += faceVolumes[i];


	std::cout << std::endl << "----------" << std::endl;
	std::cout << "the total volume of the quad mesh is " << totalVolume << std::endl;
	std::cout << std::endl << "----------" << std::endl;


}

template<size_t _N>
void tileQuad (const po::variables_map &args)
{
	
	// read in the quadMesh
    PolyMesh pmesh;
    readQuadMesh(args["quad"].as<string>(), pmesh);


    scaleQuadMesh(args["scale"].as<double>(), pmesh);


	// tile the quad mesh using the parameterTable ... 
    auto wi = WireInflator2D::construct(args["pattern"].as<string>(), args["sym"].as<int>());
    WireInflator2D::OutMeshType mesh;

    CellParameters p_params = wi->createParameters();
 
 	std::vector<CellParameters> quadParams;
    size_t nParams = p_params.numberOfParameters();
	
	vector<vector<Real>> parameterTable;
	readTable(args["inParam"].as<string>(), parameterTable);

	cout << pmesh.face.size() << "\t" << parameterTable.size() << "\t" << parameterTable[0].size() << endl;


	if (parameterTable.size() == 1)
		for(size_t i = 0; i < pmesh.face.size() - 1; ++i)
			parameterTable.push_back(parameterTable[0]);

	// check to make sure you have correct number of parameters
	if (parameterTable.size()!= pmesh.face.size()){
		throw("invalid number of rows in the parameter table!");
		return;
	} 
	else if (parameterTable[0].size() != nParams){
		throw("incorrect number of parameters!");
		return;
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
	
	// writing the final tiled mesh and adding the necessary fields to each element
	MSHFieldWriter writer(tiledMeshName, outVertices, outElements,  MeshIO::MESH_GUESS, false);


	std::vector<size_t> ids; // holds the quad id that each element belongs to
	std::vector<Real>   vols;// holds the area of the quad that each elements belongs to
	findQuadIdsForElements(outVertices, outElements, pmesh, ids, vols);


	ScalarField<Real> idsField(ids.size());
	ScalarField<Real> volsField(vols.size());
	for (size_t i = 0; i < ids.size(); ++i){
		idsField[i]  = ids[i] * 1.0;
		volsField[i] = vols[i];
	}


	// adding the id and vol elements
    writer.addField("id", idsField, DomainType::PER_ELEMENT);
    writer.addField("vol", volsField, DomainType::PER_ELEMENT);

    // adding the optimization costs 
    if (args["inCost"].as<string>() != ""){
	    vector<vector<Real>> optimizationCosts;
	    readTable(args["inCost"].as<string>(), optimizationCosts);
	    if (optimizationCosts.size() != pmesh.face.size())
	    	throw ("invalid number of rows in the cost table!");
	   
	   	ScalarField<Real> initCost(outElements.size());
	   	ScalarField<Real> finaCost(outElements.size());

	   	for (size_t i = 0; i < outElements.size(); ++i)
		{
			initCost[i] = optimizationCosts[ids[i]][0];
			finaCost[i] = optimizationCosts[ids[i]][1];
		}

		writer.addField("initCost", initCost, DomainType::PER_ELEMENT);
		writer.addField("finaCost", finaCost, DomainType::PER_ELEMENT);
	}
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
