#include "MidplaneMesher.hh"

#include <vector>
#include <utility>

#include "WireMesh.hh"
#include "PatternSignedDistance.hh"
#include "MarchingSquares/MarchingSquaresStitch.hh"
#include <filters/CurveCleanup.hh>
#include <Triangulate.h>
#include <PeriodicBoundaryMatcher.hh>
#include <stdexcept>

using namespace std;

template<class VolumeSDF>
class MidplaneSlice {
public:
    typedef typename VolumeSDF::Real   Real;
    MidplaneSlice(const VolumeSDF &vsdf)
        : m_volumeSDF(vsdf)
    {
        auto bb = vsdf.boundingBox();
        m_2DBBox = BBox<Point2<Real>>(Point2<Real>(bb.minCorner[0], bb.minCorner[1]),
                                      Point2<Real>(bb.maxCorner[0], bb.maxCorner[1]));
    }

    Point3<Real> volumePoint(const Point2<Real> &planePoint) const {
        return Point3<Real>(planePoint[0], planePoint[1], 0.0);
    }

    const BBox<Point2<Real>> &boundingBox()    const { return m_2DBBox; }
    Real signedDistance(const Point2<Real> &p) const { return m_volumeSDF.signedDistance(volumePoint(p)); }
    bool isInside(const Point2<Real> &p)       const { return signedDistance(p) <= 0; }

private:
    const VolumeSDF &m_volumeSDF;
    BBox<Point2<Real>> m_2DBBox;
};


// Steps:
//   1) Mesh Boundary (MS grid based on sqrt(max area) target)
//   2) Re-mesh boundary based on short edge criteria
//   3) Mesh with Triangle, allowing Steiner point insertion on boundary.
//      (Requires finding holes: can guess centroid and validate with SDF.
//       Another alternative: Don't mark holes, then remove afterward using
//       SDF, but this prevents Triangle from running certain operations.)
//
//   4) Post-process:
//          Snap + reflect to full period cell.
//          Optional: Run triangle *again* with boundary fixed to improve
//          meshing around symmetry plane.
template<class SignedDistanceFunction>
void MidplaneMesher<SignedDistanceFunction>::
mesh(const SignedDistanceFunction &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &triangles)
{
    auto bb = sdf.boundingBox();
    MarchingSquaresGrid msquares(meshingOptions.msGridSizeFromMaxArea(bb.dimensions()[0]),
                                 meshingOptions.msGridSizeFromMaxArea(bb.dimensions()[1]));

    MidplaneSlice<SignedDistanceFunction> slice(sdf);

#if 0
    {
        msquares.outputSignedDistanceField("sdf.msh", slice);
    }
#endif

    // Get ccw-ordered segments of boundary/interior edges
    auto result = msquares.extractBoundaryPolygons(slice, 0.0);

    std::vector<std::pair<size_t, size_t>> edges;
    edges.reserve(result.numEdges());
    for (const auto &s : result.segments)
        edges.insert(edges.end(), s.second.begin(), s.second.end());

    Real maxLen = meshingOptions.maxEdgeLenFromMaxArea();
    Real minLen = meshingOptions.minEdgeLenFromMaxArea(std::min(bb.dimensions()[0],
                                                                bb.dimensions()[1]));
    // std::cout << "Using maxLen " << maxLen << ", minLen " << minLen << std::endl;

    // Organized polygon soup into ccw polygons
    std::list<std::list<Point2D>> polygons;
    extract_polygons<2>(result.points, edges, polygons);

    // Non-periodic polygon cleanup since we're meshing a quarter-cell.
    // Shouldn't matter anyway because the periodic should only have long edges.
    for (auto &poly : polygons) {
        curveCleanup<2>(poly, slice.boundingBox(), minLen, maxLen,
                meshingOptions.featureAngleThreshold, false);
    }

#if 0
    {
        IOElementEdgeSoupFromClosedPolygonList<Point2D> esoup(polygons);
        MeshIO::save("cleaned_polygons.msh", esoup);
    }
#endif

    // Determine which polygon is touching the bbox (there must be exactly one):
    // this is the only non-hole polygon.
    std::vector<bool> isHoleBdry;
    size_t numHoles = 0;
    for (const auto &poly : polygons) {
        bool isHole = true;
        for (const auto &p : poly) {
            if (PeriodicBoundaryMatcher::FaceMembership<2>(p, slice.boundingBox()).count()) {
                isHole = false;
                break;
            }
        }
        isHoleBdry.push_back(isHole);
        numHoles += isHole;
    }
    if (polygons.size() - numHoles != 1) {
        throw std::runtime_error("Should have exactly one bbox-incident curve; got "
                + std::to_string(polygons.size() - numHoles) + ".");
    }

    // Try to find a point inside each hole boundary by offsetting from the
    // boundary curve and checking if the point is outside the object.
    std::vector<Point2D> holePts;
    {
        size_t i = 0;
        for (const auto &poly : polygons) {
            if (poly.size() < 3) throw std::runtime_error("Polygon of size " + std::to_string(poly.size()) + " in marching squares output.");
            if (isHoleBdry.at(i++)) {
                bool found = false;
                for (auto p_it = poly.begin(); !found && (p_it != poly.end()); ++p_it) {
                    auto p_next = p_it; ++p_next;
                    if (p_next == poly.end()) p_next = poly.begin();
                    Point2D midpoint = 0.5 * (*p_next + *p_it);
                    Point2D edgeVec = *p_next - *p_it; // ccw pointing edge
                    Vector2D inwardNormal(-edgeVec[1], edgeVec[0]);
                    Point2D query = midpoint + (1e-5 / inwardNormal.norm()) * inwardNormal;
                    if (slice.signedDistance(query) < 0.0) {
                        holePts.push_back(query);
                        found = true;
                    }
                }
                if (!found) {
                    throw std::runtime_error("Couldn't find point inside hole " + std::to_string(i));
                }
            }
        }
    }
    if (polygons.size() - holePts.size() != 1) {
        throw std::runtime_error("Should have exactly one non-hole curve. Got "
                + std::to_string(polygons.size() - holePts.size()));
    }

    triangulatePSLC(polygons, holePts, vertices, triangles,
                    meshingOptions.maxArea, "Q");
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
template class MidplaneMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>>>;
