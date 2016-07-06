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

#define     DEBUG_OUT 0
#define SDF_DEBUG_OUT 0

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
//  Potential future improvement: mesh full cell by reflecting and remeshing the
//  interior.
template<class SignedDistanceFunction>
void MidplaneMesher<SignedDistanceFunction>::
mesh(const SignedDistanceFunction &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &triangles)
{
    auto bb = sdf.boundingBox();
    MarchingSquaresGrid msquares(meshingOptions.msGridSizeFromMaxArea(bb.dimensions()[0]),
                                 meshingOptions.msGridSizeFromMaxArea(bb.dimensions()[1]),
                                 meshingOptions.marchingSquaresCoarsening);

    MidplaneSlice<SignedDistanceFunction> slice(sdf);

#if SDF_DEBUG_OUT
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

    Real maxLen = meshingOptions.maxBdryEdgeLen();
    Real domainLength = std::min(bb.dimensions()[0], bb.dimensions()[1]);
    Real minLen = meshingOptions.minEdgeLenFromMaxArea(domainLength);
    // std::cout << "Using maxLen " << maxLen << ", minLen " << minLen << "old maxLen: " << meshingOptions.maxEdgeLenFromMaxArea() << std::endl;

    // Organized polygon soup into ccw polygons
    std::list<std::list<Point2D>> polygons;
    extract_polygons<2>(result.points, edges, polygons);

    // Non-periodic polygon cleanup (we assume we're meshing a quarter-cell).
    // TODO: periodic cleanup if we're actually meshing the full period cell.
    std::vector<std::vector<Real>> variableMinLens;
    variableMinLens.reserve(polygons.size());
    bool periodic = false;
    if (meshingOptions.curvatureAdaptive) {
        for (const auto &poly : polygons) {
            std::vector<MeshIO::IOVertex> vertices;
            std::vector<MeshIO::IOElement> elements;
            std::vector<Real> signedCurvatures;

            size_t offset = vertices.size();
            vertices.reserve(offset + poly.size());
            for (const auto &p : poly) {
                vertices.emplace_back(p);
                elements.emplace_back(vertices.size() - 1, vertices.size());
            }
            elements.back()[1] = offset;
            signedCurvatures.reserve(vertices.size());

            // Use a finite difference to determine if the normal points
            // left or right from a curve tangent. But make sure we test a true
            // boundary segment (instead of a cell boundary segment)
            auto segmentP0 = poly.begin();
            for (; segmentP0 != poly.end(); ++segmentP0) {
                if (!PeriodicBoundaryMatcher::FaceMembership<2>(*segmentP0, slice.boundingBox()).onAnyFace())
                    break;
            }

            // Skip all-cell-boundary (or an annoying corner case of nearly
            // all-cell-boundary) polygons--these don't need adaptive meshing
            auto segmentP1 = segmentP0;
            if ((segmentP0 == poly.end()) || (++segmentP1 == poly.end())) {
                variableMinLens.emplace_back();
                continue;
            }

            // Choose curvature sign so that it reflects concave/convex geometry
            // (object normal is a 90 clockwise rotation of curve tangent)
            Point2D midpoint = 0.5 * (*segmentP0 + *segmentP1);
            Vector2D tangent = *segmentP1 - *segmentP0;
            tangent *= 1.0 / tangent.norm();
            Vector2D right(tangent[1], -tangent[0]);
            Real eps = 1e-3;
            // Positive if "right" is an outward normal
            Real diff = slice.signedDistance(midpoint + eps * right) -
                        slice.signedDistance(midpoint - eps * right);
            Real sign = diff > 0;

            // With this sign convention, convex geometry
            // (tangent turning ccw towards interior) has positive curvature and
            // concave geometry (tangent turning cw towards exterior) has
            // negative sign.
            for (Real k : signedCurvature(poly))
                signedCurvatures.push_back(sign * k);

            assert(signedCurvatures.size() == vertices.size());
            ScalarField<Real> kappa(vertices.size());
            for (size_t i = 0; i < vertices.size(); ++i)
                kappa[i] = signedCurvatures[i];

#if DEBUG_OUT
            MSHFieldWriter writer("curvatures.msh", vertices, elements);
            writer.addField("signed curvature", kappa, DomainType::PER_NODE);
#endif
            // {
            //     std::ofstream curvatureOut("curvature.txt");
            //     for (size_t i = 0; i < signedCurvatures.size(); ++i)
            //         curvatureOut << signedCurvatures[i] << std::endl;
            // }

            // Chose adaptive edge length based on curvature: highly negative
            // curvature uses the fine marching squares-based edge length while the
            // zero and higher curvature uses maxLen / 4
            ScalarField<Real> lengths(vertices.size()); // Actually per edge
            for (size_t i = 0; i < vertices.size(); ++i) {
                auto pt1 = truncateFrom3D<Point2D>(vertices[i]),
                     pt2 = truncateFrom3D<Point2D>(vertices[(i + 1) % vertices.size()]);
                int numAverage = 0;
                Real k = 0;
                if (!PeriodicBoundaryMatcher::FaceMembership<2>(pt1, slice.boundingBox()).onAnyFace()) {
                    ++numAverage; k += kappa[i];
                }
                if (!PeriodicBoundaryMatcher::FaceMembership<2>(pt2, slice.boundingBox()).onAnyFace()) {
                    ++numAverage; k += kappa[(i + 1) % vertices.size()];
                }
                if (numAverage != 0) k /= numAverage;

                // Want to interpolate from upper at k >= c to lower at k = d
                // using function a * 2^(k * b)
                Real upper = maxLen / 4.0;
                Real lower = minLen;
                Real c = -1.0, d = -4;
                // upper = a * 2^(cb)
                // lower = a * 2^(db)
                // upper / lower = 2^((c - d) b) ==> b = log2(u / l) / (c - d)
                // a = upper / 2^(c * b)
                Real b = log(upper / lower) / (log(2) * (c - d));
                Real a = upper / pow(2, c * b);
                lengths[i] = a * pow(2, b * k);
                lengths[i] = std::min(lengths[i], upper);
                lengths[i] = std::max(lengths[i], lower);
            }
#if DEBUG_OUT
            writer.addField("min_lengths", lengths, DomainType::PER_ELEMENT);
            ScalarField<Real> edgeLengths(vertices.size());
            ScalarField<Real> isShort(vertices.size());
            for (size_t i = 0; i < vertices.size(); ++i) {
                auto pt1 = truncateFrom3D<Point2D>(vertices[i]),
                     pt2 = truncateFrom3D<Point2D>(vertices[(i + 1) % vertices.size()]);
                edgeLengths[i] = (pt2 - pt1).norm();
                isShort[i] = (edgeLengths[i] < lengths[i]) ? 1.0 : 0.0;
            }

            writer.addField("edge_lengths", edgeLengths, DomainType::PER_ELEMENT);
            writer.addField("is_short", isShort, DomainType::PER_ELEMENT);
#endif
            variableMinLens.emplace_back(lengths.domainSize());
            std::vector<Real> &vml = variableMinLens.back();
            for (size_t i = 0; i < lengths.domainSize(); ++i)
                vml[i] = lengths[i];
        }
    }

#if DEBUG_OUT
    std::cout << polygons.size() << " polygons. Sizes:" << std::endl;
    for (auto &poly : polygons) {
        std::cout << "\t" << poly.size() << std::endl;
    }

    {
        IOElementEdgeSoupFromClosedPolygonList<Point2D> esoup(polygons);
        MeshIO::save("ms_polygons.msh", esoup);
    }
#endif

    BENCHMARK_START_TIMER("Curve Cleanup");
    {
        size_t i = 0;
        for (auto &poly : polygons) {
            Real cellEpsilon = 0.0; // Marching squares guarantees cell boundary vertex coords are exact
            if (variableMinLens.size()) {
                curveCleanup<2>(poly, slice.boundingBox(), minLen, maxLen,
                        meshingOptions.featureAngleThreshold, periodic, variableMinLens.at(i), cellEpsilon);
            }
            else {
                curveCleanup<2>(poly, slice.boundingBox(), minLen, maxLen,
                        meshingOptions.featureAngleThreshold, periodic, std::vector<Real>(), cellEpsilon);
            }
            ++i;
        }
    }
    BENCHMARK_STOP_TIMER("Curve Cleanup");

#if DEBUG_OUT
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
                // Choose the hole point of greatest signed distance (furthest
                // into void) for robustness
                Real maxDist = 0;
                Point2D candidate;
                for (auto p_it = poly.begin(); p_it != poly.end(); ++p_it) {
                    auto p_next = p_it; ++p_next;
                    if (p_next == poly.end()) p_next = poly.begin();
                    Point2D midpoint = 0.5 * (*p_next + *p_it);
                    Point2D edgeVec = *p_next - *p_it;
                    Vector2D outwardNormal(edgeVec[1], -edgeVec[0]); // outward from geometry (inward to hole)
                    Point2D query = midpoint + (1e-4 / outwardNormal.norm()) * outwardNormal;
                    Real querySD = slice.signedDistance(query);
                    if (querySD > maxDist) {
                        maxDist = querySD;
                        candidate = query;
                    }
                }
                if (maxDist == 0) throw std::runtime_error("Couldn't find point inside hole " + std::to_string(i));
                holePts.push_back(candidate);
                // std::cerr << "Found hole point: " << candidate << std::endl;
                // std::cerr << "signed distance at hole point: " << maxDist << std::endl;
            }
        }
    }
    if (polygons.size() - holePts.size() != 1) {
        throw std::runtime_error("Should have exactly one non-hole curve. Got "
                + std::to_string(polygons.size() - holePts.size()));
    }

    triangulatePSLC(polygons, holePts, vertices, triangles,
                    meshingOptions.maxArea, "Q");

#if DEBUG_OUT
    MeshIO::save("triangulated_polygon.msh", vertices, triangles);
#endif

}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
template class MidplaneMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>>>;
