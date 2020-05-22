//
// Created by Davi Colli Tozoni on 5/11/18.
//

#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/LinearElasticity.hh>

#include "shape_optimization/PolygonMesher.hh"
#include <MeshFEM/MSHFieldWriter.hh>
#include <isosurface_inflator/MeshingOptions.hh>

#include <iostream>
#include <boost/algorithm/string.hpp>

#include <boost/program_options.hpp>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/Parallelism.hh>

#include <inflators/wrappers/BoundaryPerturbationInflator.hh>

#include <MeshFEM/Utilities/EdgeSoupAdaptor.hh>

#include <shape_optimization/ParametersMask.hh>

namespace po = boost::program_options;
using namespace std;

template<size_t _N>
using HMG = LinearElasticity::HomogenousMaterialGetter<Materials::Constant>::template Getter<_N>;
typedef LinearElasticity::Mesh<2, 2, HMG> Mesh;
typedef LinearElasticity::Simulator<Mesh> Simulator;

void usage(int status, const po::options_description &visible_opts) {
    cerr << "Usage: ./polygonMesher_cli <polygons file> out.msh" << endl;
    cerr << "eg: ./polygonMesher_cli <polygons file> out.msh" << endl;
    cout << visible_opts << endl;
    exit(status);
}

po::variables_map parseCmdLine(int argc, char *argv[]) {
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
            ("input", po::value<string>(),  "input polygons/mesh file")
            ("outMSH", po::value<string>(), "output msh file")
            ;
    po::positional_options_description p;
    p.add("input", 1);
    p.add("outMSH", 1);

    // Options visible in the help message.
    po::options_description visible_opts;
    visible_opts.add_options()("help,h", "Produce this help message")
            ("mopts,m", po::value<string>(),               "Meshing options file")
            ("periodic",                                   "makes the mesh tileable")
            ("polygonInput",                               "defines that the input file is a set of polygons instead of a mesh")
            ("remeshingRegion,r", po::value<string>(),     "file describing (non)remeshing regions")
            ("gridSize,g", po::value<size_t>(),            "grid size to be used in resampling")
            ("noneCurveCleanup,n",                         "do not clean up boundary of polygon")
            ("collapseCleanup,c",                          "perform collapse clean up")
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
    if (vm.count("input") == 0) {
        cout << "Error: must specify boundary polygons" << endl;
        fail = true;
    }

    if (vm.count("outMSH") == 0) {
        cout << "Error: must specify output" << endl;
        fail = true;
    }

    if (fail)
        usage(1, visible_opts);

    return vm;
}

std::list<std::list<Point2D>> readPolygons(string polygonsPath) {
    std::vector<MeshIO::IOVertex> vertices;
    std::vector<MeshIO::IOElement> elements;
    std::list<std::list<Point2D>> polygons;

    MeshIO::load(polygonsPath, vertices, elements);

    //loop through elements (polygons)
    for (auto element : elements) {
        std::list<Point2D> poly;

        //loop through vertices of element
        for (auto vertexIdx : element) {
            poly.push_back(vertices[vertexIdx]);
        }

        polygons.push_back(poly);
    }

    return polygons;
}

std::list<std::list<Point2D>> readMesh(string meshPath, PolygonMesher &polygonMesher) {
    std::list<std::list<Point2D>> polygons;
    std::vector<MeshIO::IOVertex> meshVertices;
    std::vector<MeshIO::IOElement> meshElements;
    MeshIO::load(meshPath, meshVertices, meshElements);

    Simulator sim(meshElements, meshVertices);

    std::vector<MeshIO::IOVertex> boundaryVertices;
    std::vector<MeshIO::IOElement> boundaryElements;
    polygonMesher.meshToBoundaryMeshIO(sim.mesh(), boundaryVertices, boundaryElements);
    polygons = polygonMesher.extractPolygons(boundaryVertices, boundaryElements);

    //std::cout << "polygons: " << polygons.size() << std::endl;

    return polygons;
}

void extractFilteringRegions(string bcondsPath, std::list<std::list<Point2D>> polygons, vector<std::unique_ptr<Region<Point3D>>> &remeshingRegions, vector<std::unique_ptr<Region<Point3D>>> &exceptRegions) {
    vector<Point3D> points;

    for (auto poly : polygons) {
        for (auto vertex : poly)
            points.push_back(padTo3D(vertex));
    }

    if (!bcondsPath.empty()) {
        BBox<Point3D> bb(points);
        auto dim = bb.dimensions();
        if ((std::abs(dim[0]) < 1e-6) || (std::abs(dim[1]) < 1e-6))
            throw std::runtime_error("Degenerate pattern");

        bool no_rigid_motion;
        vector<std::unique_ptr<Region<Point3D>>> glueRegions;
        ParametersMask::jsonToRegions(bcondsPath, bb, remeshingRegions, exceptRegions, glueRegions);

        // In remeshing case, glue regions should also be not remeshed, so append all its terms to exceptRegions
        // Do it in a manual way, because of unique_ptr
        for (size_t i = 0; i < glueRegions.size(); i++) {
            exceptRegions.emplace_back(std::move(glueRegions[i]));
        }
        //exceptRegions.insert(exceptRegions.end(), glueRegions.begin(), glueRegions.end());
    }
}

int main(int argc, char *argv[])
{
    std::list<std::list<Point2D>> polygons;
    std::vector<MeshIO::IOVertex> vertices;
    std::vector<MeshIO::IOElement> elements;

    auto args = parseCmdLine(argc, argv);

    PolygonMesher polygonMesher;

    if (args.count("mopts"))
        polygonMesher.meshingOptions.load(args["mopts"].as<string>());

    if (args.count("polygonInput") > 0)
        polygons = readPolygons(args["input"].as<string>());
    else
        polygons = readMesh(args["input"].as<string>(), polygonMesher);

    vector<std::unique_ptr<Region<Point3D>>> remeshingRegions;
    vector<std::unique_ptr<Region<Point3D>>> exceptRegions;
    if (args.count("remeshingRegion") > 0) {
        extractFilteringRegions(args["remeshingRegion"].as<string>(), polygons, remeshingRegions, exceptRegions);
    }

    if (args.count("gridSize") > 0) {
        polygonMesher.meshingOptions.marchingSquaresGridSize = args["gridSize"].as<size_t>();
    }

    if (args.count("periodic") > 0) {
        polygonMesher.meshInterfaceConsistently = true;
    }

    if (args.count("noneCurveCleanup") > 0) {
        polygonMesher.meshingOptions.curveSimplifier = MeshingOptions::NONE;
    }
    else if (args.count("collapseCleanup") > 0) {
        polygonMesher.meshingOptions.curveSimplifier = MeshingOptions::COLLAPSE;
    }

    polygonMesher.mesh(polygons, vertices, elements, remeshingRegions, exceptRegions);

    MeshIO::save(args["outMSH"].as<string>(), vertices, elements);

    return 0;
}
