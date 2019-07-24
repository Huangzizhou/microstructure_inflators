//
// Created by Davi Colli Tozoni on 5/11/18.
//

#include "PolygonMesher.hh"

#include <vector>
#include <utility>

#include <MeshFEM/filters/CurveCleanup.hh>
#include <MeshFEM/filters/ResampleCurve.hh>
#include <MeshFEM/Triangulate.h>
#include <MeshFEM/PeriodicBoundaryMatcher.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <stdexcept>

#include <MeshFEM/Utilities/EdgeSoupAdaptor.hh>

#include <nonstd/optional.hpp>

#define     DEBUG_OUT 1

std::list<std::list<Point2D>> PolygonMesher::
extractPolygons(std::vector<MeshIO::IOVertex> &vertices, std::vector<MeshIO::IOElement> &elements) const {
    std::list<std::list<Point2D>> polygons;
    std::vector<bool> used(vertices.size(), false);
    std::vector<std::vector<int>> neighbor(vertices.size());

    //loop through elements
    for (auto element : elements) {
        neighbor[element[0]].push_back(element[1]);
        neighbor[element[1]].push_back(element[0]);
    }

    for (size_t initialVertex = 0; initialVertex < vertices.size(); initialVertex++) {
        if (used[initialVertex])
            continue;

        std::list<Point2D> newPolygon;
        size_t v = initialVertex;
        while (!used[v]) {
            newPolygon.push_back(vertices[v]);
            used[v] = true;

            assert(neighbor[v].size() == 2);

            v = used[neighbor[v][0]] ? neighbor[v][1] : neighbor[v][0];
        }

        polygons.push_back(newPolygon);
    }

    return polygons;
}

void PolygonMesher::
mesh(std::list<std::list<Point2D>> &polygons,
     std::vector<MeshIO::IOVertex>  &vertices,
     std::vector<MeshIO::IOElement> &triangles,
     const std::vector<std::unique_ptr<Region<Point3D>>> &remeshingRegions,
     const std::vector<std::unique_ptr<Region<Point3D>>> &exceptRegions) const
{
    std::list<Point2D> outerBoundary = polygons.front();
    BBox<Point2d> bb(outerBoundary);
    size_t gridSizeX = meshingOptions.msGridSizeFromMaxArea(bb.dimensions()[0]);
    double gridCellWidth = bb.dimensions()[0] / gridSizeX;
    double maxLen = meshingOptions.maxBdryEdgeLen();
    double domainLength = std::min(bb.dimensions()[0], bb.dimensions()[1]);
    double minLen = meshingOptions.minEdgeLenFromMaxArea(domainLength);
    // std::cout << "Using maxLen " << maxLen << ", minLen " << minLen << "old maxLen: " << meshingOptions.maxEdgeLenFromMaxArea() << std::endl;

    // Remove extremely small polygons (relative to MS grid) that are probably
    // due to sampling error, and that will cause robustness issues for hole
    // finding.
    polygons.remove_if([gridCellWidth](const std::list<Point2D> &poly) {
                           BBox<Point2D> pbb(poly);
                           auto dims = pbb.dimensions();
                           return ((dims[0] < 0.25 * gridCellWidth) &&
                                   (dims[1] < 0.25 * gridCellWidth));
                       }
    );

    using PolygonAdaptor = IOElementEdgeSoupFromClosedPolygonCollection<decltype(polygons)>;

#if DEBUG_OUT
    std::cout << polygons.size() << " polygons. Sizes:" << std::endl;
    for (auto &poly : polygons)
        std::cout << "\t" << poly.size() << std::endl;

    MeshIO::save("ms_polygons.msh", PolygonAdaptor(polygons));
#endif

    if (!msPolygonPath.empty())
        MeshIO::save(msPolygonPath, PolygonAdaptor(polygons));

    BENCHMARK_START_TIMER("Curve Cleanup");
    {
        // Only remesh the cell boundary if we need to "meshInterfaceConsistently"
        // (i.e. prevent Triangle from inserting Steiner points). Otherwise,
        // we'll let Triangle subdivide the boundary optimally.
        nonstd::optional<double> cellBdryEdgeLen;
        if (this->meshInterfaceConsistently) cellBdryEdgeLen = meshingOptions.maxEdgeLenFromMaxArea() / 2;

        if (meshingOptions.curveSimplifier == MeshingOptions::COLLAPSE) {
            size_t i = 0;
            for (auto &poly : polygons) {
                curveCleanup<2>(poly, bb, minLen, maxLen,
                                meshingOptions.featureAngleThreshold, this->meshInterfaceConsistently, cellBdryEdgeLen,
                                std::vector<double>(), 0.0);
                ++i;
            }
        }
        else if (meshingOptions.curveSimplifier == MeshingOptions::RESAMPLE) {
            size_t i = 0;
            Real targetLen = minLen * 3;
            for (auto &poly : polygons) {
                resampleCurve<2>(poly, bb, targetLen,
                                 meshingOptions.featureAngleThreshold,
                                 remeshingRegions,
                                 exceptRegions,
                                 cellBdryEdgeLen,
                                 std::vector<double>(),
                                 0.0, !this->meshInterfaceConsistently);
                ++i;
            }
        }
        else if (meshingOptions.curveSimplifier == MeshingOptions::NONE) {
            // does not do anything
        }
        else throw std::runtime_error("Illegal curve simplifier");
    }
    BENCHMARK_STOP_TIMER("Curve Cleanup");

    // Remove all invalid polygons
    for (auto it = polygons.begin(); it != polygons.end(); ) {
        if (it->size() > 2) {
            ++it;
        } else {
            it = polygons.erase(it);
        }
    }

#if DEBUG_OUT
    std::cout << "Simplified polygon sizes:" << std::endl;
    for (auto &poly : polygons)
        std::cout << "\t" << poly.size() << std::endl;
    MeshIO::save("cleaned_polygons.msh",
                 PolygonAdaptor(polygons));
#endif

    // Determine which polygon is touching the bbox (there should be exactly one):
    // this is the only non-hole polygon.
    // Using the actuall bounding box of polygons vertices should take care of cases
    // including non periodic instances
    std::vector<bool> isHoleBdry;

    // Find bounding box of entire shape
    Real leftist = std::numeric_limits<Real>::max();
    BBox<Point2d> cleanedbb;
    for (const auto &poly : polygons) {
        BBox<Point2d> candidate(poly);
        if (candidate.minCorner[0] < leftist) {
            leftist = candidate.minCorner[0];
            cleanedbb = candidate;
        }

    }
    //BBox<Point2d> cleanedbb(polygons.front());
    size_t numHoles = 0;
    double tol = 1e-5;
    for (const auto &poly : polygons) {
        bool isHole = true;
        for (const auto &p : poly) {
            if (PeriodicBoundaryMatcher::FaceMembership<2>(p, cleanedbb, tol).count()) {
                isHole = false;
                break;
            }
        }
        isHoleBdry.push_back(isHole);
        numHoles += isHole;
    }

    // Actually, we can have more than one bbox-incident curve when a
    // neigboring cell's geometry extends into this cell.
#if 0
    if (polygons.size() - numHoles != 1) {
        throw std::runtime_error("Should have exactly one bbox-incident curve; got "
                + std::to_string(polygons.size() - numHoles) + ".");
    }
#endif

    // Try to find a point inside each hole boundary.
    BENCHMARK_START_TIMER("Hole detection");
    std::vector<Point2D> holePts;
    {
        size_t i = 0;
        for (const auto &poly : polygons) {

            if (poly.size() < 3) throw std::runtime_error("Polygon of size " + std::to_string(poly.size()) + " in marching squares output.");
            if (isHoleBdry.at(i++)) {
                // Brute-force solution to robustly finding point in the hole:
                // Triangulate hole and then consider triangle barycenters

                //TODO: possible to improve performance by an algorithm that computes directly a point in the interior
                // of the hole

                std::list<std::list<Point2D>> holePolygons(1, poly);
                std::vector<MeshIO::IOVertex > holeVertices;
                std::vector<MeshIO::IOElement> holeTriangles;
                //MeshIO::save("polygons.msh", PolygonAdaptor(polygons));
                triangulatePSLC(holePolygons, std::vector<Point2D>(),
                                holeVertices, holeTriangles, 1.0, "Q");

                // Choose barycenter of first triangle
                MeshIO::IOElement tri = holeTriangles[0];
                Point2D holePoint = truncateFrom3D<Point2D>(
                        1.0 / 3.0 * (holeVertices[tri[0]].point +
                                     holeVertices[tri[1]].point +
                                     holeVertices[tri[2]].point));

                holePts.push_back(holePoint);
#if DEBUG_OUT
                std::cerr << "Found hole point: " << holePoint << std::endl;
#endif
            }
        }
    }

    if (holePts.size() != numHoles) {
        std::cerr << "WARNING: couldn't find all holes" << std::endl;
        std::cerr << "WARNING: There were suppose to be " << numHoles
                  << " holes. We found only: " << holePts.size() << std::endl;
    }
    BENCHMARK_STOP_TIMER("Hole detection");


    MeshIO::save("polygons.msh", PolygonAdaptor(polygons));

    if (polygons.size() == 0) return;
    triangulatePSLC(polygons, holePts, vertices, triangles,
                    meshingOptions.maxArea,
                    (this->meshInterfaceConsistently ? "QY" : "Q"));

#if DEBUG_OUT
    std::cout << "# vertices: " << vertices.size() << std::endl;
    std::cout << "# triangles: " << triangles.size() << std::endl;
    MeshIO::save("triangulated_polygon.msh", vertices, triangles);
#endif
}
