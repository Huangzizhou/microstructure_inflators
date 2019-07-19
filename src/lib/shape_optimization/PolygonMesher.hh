//
// Created by Davi Colli Tozoni on 5/11/18.
//

#ifndef MICROSTRUCTURES_POLYGONMESHER_H
#define MICROSTRUCTURES_POLYGONMESHER_H

#include <MeshFEM/MeshIO.hh>
#include <isosurface_inflator/MeshingOptions.hh>
#include <list>
#include <string>

class PolygonMesher {
public:
    PolygonMesher(bool keepBoundaryIntact = true) {
        meshInterfaceConsistently = keepBoundaryIntact;
    }

    void mesh(std::list<std::list<Point2D>> &polygons,
              std::vector<MeshIO::IOVertex> &vertices,
              std::vector<MeshIO::IOElement> &elements,
              const std::vector<std::unique_ptr<Region<Point3D>>> &remeshingRegions = std::vector<std::unique_ptr<Region<Point3D>>>(),
              const std::vector<std::unique_ptr<Region<Point3D>>> &exceptRegions = std::vector<std::unique_ptr<Region<Point3D>>>()) const;

    template<class _Mesh>
    void meshToBoundaryMeshIO(const _Mesh &mesh, std::vector<MeshIO::IOVertex> &outVertices, std::vector<MeshIO::IOElement> &outElements) const {
        outElements.resize(mesh.numBoundaryElements());
        for (size_t bei = 0; bei < mesh.numBoundaryElements(); ++bei) {
            auto be = mesh.boundaryElement(bei);
            for (size_t c = 0; c < be.numVertices(); ++c) {
                outElements[bei].push_back(be.vertex(c).index());
            }
        }

        outVertices.reserve(mesh.numBoundaryVertices());
        for (size_t bvi = 0; bvi < mesh.numBoundaryVertices(); ++bvi) {
            Point2D point = mesh.boundaryVertex(bvi).node().volumeNode()->p;
            outVertices.push_back(point);
        }
    }

    std::list<std::list<Point2D>> extractPolygons(std::vector<MeshIO::IOVertex> &vertices,
                                                  std::vector<MeshIO::IOElement> &elements) const;

    std::string msPolygonPath;
    MeshingOptions meshingOptions;
    bool meshInterfaceConsistently;
};

#endif //MICROSTRUCTURES_POLYGONMESHER_H
