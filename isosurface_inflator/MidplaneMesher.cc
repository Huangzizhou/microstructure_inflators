#include "MidplaneMesher.hh"

#include <vector>
#include <utility>

#include "MarchingSquares/MarchingSquaresStitch.hh"
#include <filters/CurveCleanup.hh>
#include <Triangulate.h>
#include <PeriodicBoundaryMatcher.hh>
#include <stdexcept>

#define     DEBUG_OUT 0
#define SDF_DEBUG_OUT 0

using namespace std;

struct MidplaneSlice : public SignedDistanceRegion<2> {
    MidplaneSlice(const SignedDistanceRegion<3> &volumeRegion)
        : m_volRegion(volumeRegion) {
        auto bb = volumeRegion.boundingBox();
        m_2DBBox = BBox<Point2d>(Point2d(bb.minCorner[0], bb.minCorner[1]),
                                 Point2d(bb.maxCorner[0], bb.maxCorner[1]));
    }

    Point3d volumePoint(const Point2d &planePoint) const {
        return Point3d(planePoint[0], planePoint[1], 0.0);
    }

    virtual const BBox<Point2d> &boundingBox()      const override { return m_2DBBox; }
    virtual double signedDistance(const Point2d &p) const override { return m_volRegion.signedDistance(volumePoint(p)); }

private:
    const SignedDistanceRegion<3> &m_volRegion;
    BBox<Point2d> m_2DBBox;
};

void computeCurvatureAdaptiveMinLen(const std::list<std::list<Point2D>> &polygons,
                                    const MidplaneSlice &slice,
                                    const double minLen, const double maxLen,
                                    std::vector<std::vector<double>> &variableMinLens);

// Steps:
//   1) Mesh Boundary (MS grid based on sqrt(max area) target)
//   2) Re-mesh boundary based on short edge criteria
//   3) Mesh with Triangle, allowing Steiner point insertion on boundary.
//      (Requires finding holes: can guess centroid and validate with SDF.
//       Another alternative: Don't mark holes, then remove afterward using
//       SDF, but this prevents Triangle from running certain operations.)
//  Potential future improvement: mesh full cell by reflecting and remeshing the
//  interior.
void MidplaneMesher::
mesh(const SignedDistanceRegion<3> &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &triangles)
{
    auto bb = sdf.boundingBox();
    size_t gridSizeX = meshingOptions.msGridSizeFromMaxArea(bb.dimensions()[0]),
           gridSizeY = meshingOptions.msGridSizeFromMaxArea(bb.dimensions()[1]);
    double gridCellWidth = bb.dimensions()[0] / gridSizeX;
    MarchingSquaresGrid msquares(gridSizeX, gridSizeY,
                                 meshingOptions.marchingSquaresCoarsening);

    MidplaneSlice slice(sdf);

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

    double maxLen = meshingOptions.maxBdryEdgeLen();
    double domainLength = std::min(bb.dimensions()[0], bb.dimensions()[1]);
    double minLen = meshingOptions.minEdgeLenFromMaxArea(domainLength);
    // std::cout << "Using maxLen " << maxLen << ", minLen " << minLen << "old maxLen: " << meshingOptions.maxEdgeLenFromMaxArea() << std::endl;

    // Organized polygon soup into ccw polygons
    std::list<std::list<Point2D>> polygons;
    extract_polygons<2>(result.points, edges, polygons);

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

    // Clean up/simplify marching squares polygons.
    std::vector<std::vector<double>> variableMinLens;
    if (meshingOptions.curvatureAdaptive)
        computeCurvatureAdaptiveMinLen(polygons, slice, minLen, maxLen, variableMinLens);

#if DEBUG_OUT
    std::cout << polygons.size() << " polygons. Sizes:" << std::endl;
    for (auto &poly : polygons)
        std::cout << "\t" << poly.size() << std::endl;

    MeshIO::save("ms_polygons.msh",
                 IOElementEdgeSoupFromClosedPolygonList<Point2D>(polygons));
#endif

    BENCHMARK_START_TIMER("Curve Cleanup");
    {
        size_t i = 0;
        for (auto &poly : polygons) {
            curveCleanup<2>(poly, slice.boundingBox(), minLen, maxLen,
                    meshingOptions.featureAngleThreshold, this->periodic,
                    variableMinLens.size() ? variableMinLens.at(i) : std::vector<double>(),
                    0.0 /* marching squares guarantees cell boundary vertex coords are exact */);
            ++i;
        }
    }
    BENCHMARK_STOP_TIMER("Curve Cleanup");

#if DEBUG_OUT
    MeshIO::save("cleaned_polygons.msh",
                 IOElementEdgeSoupFromClosedPolygonList<Point2D>(polygons));
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

                std::list<std::list<Point2D>> holePolygons(1, poly);
                std::vector<MeshIO::IOVertex > holeVertices;
                std::vector<MeshIO::IOElement> holeTriangles;
                triangulatePSLC(holePolygons, std::vector<Point2D>(),
                                holeVertices, holeTriangles, 1.0, "Q");

                double maxDist = std::numeric_limits<double>::lowest();
                Point2D bestCandidate;
                // Choose barycenter with greatest signed distance (furthest
                // into hole) for robustness
                for (const auto &tri : holeTriangles) {
                    auto candidate = truncateFrom3D<Point2D>(
                            1.0 / 3.0 * (holeVertices[tri[0]].point +
                                         holeVertices[tri[1]].point +
                                         holeVertices[tri[2]].point));
                    double sd = slice.signedDistance(candidate);
                    if (sd > maxDist) {
                        bestCandidate = candidate;
                        maxDist = sd;
                    }
                }

                if (maxDist == std::numeric_limits<double>::lowest())
                    throw std::runtime_error("Couldn't find point inside hole " + std::to_string(i));
                if (maxDist <= 0)
                    std::cerr << "WARNING: hole point inside object (nonpostive signed distance)" << std::endl;
                holePts.push_back(bestCandidate);
#if DEBUG_OUT
                std::cerr << "Found hole point: " << bestCandidate << std::endl;
                std::cerr << "signed distance at hole point: " << maxDist << std::endl;
#endif
            }
        }
    }
    BENCHMARK_STOP_TIMER("Hole detection");
    if (polygons.size() - holePts.size() != 1) {
        throw std::runtime_error("Should have exactly one non-hole curve. Got "
                + std::to_string(polygons.size() - holePts.size()));
    }

    triangulatePSLC(polygons, holePts, vertices, triangles,
                    meshingOptions.maxArea,
                    (this->periodic ? "QY" : "Q"));

#if DEBUG_OUT
    MeshIO::save("triangulated_polygon.msh", vertices, triangles);
#endif
}

void computeCurvatureAdaptiveMinLen(const std::list<std::list<Point2D>> &polygons,
                                    const MidplaneSlice &slice,
                                    const double minLen, const double maxLen,
                                    std::vector<std::vector<double>> &variableMinLens)
{
    variableMinLens.reserve(polygons.size());
    for (const auto &poly : polygons) {
        std::vector<MeshIO::IOVertex> vertices;
        std::vector<MeshIO::IOElement> elements;
        std::vector<double> signedCurvatures;

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
        double eps = 1e-3;
        // Positive if "right" is an outward normal
        double diff = slice.signedDistance(midpoint + eps * right) -
                      slice.signedDistance(midpoint - eps * right);
        double sign = diff > 0;

        // With this sign convention, convex geometry
        // (tangent turning ccw towards interior) has positive curvature and
        // concave geometry (tangent turning cw towards exterior) has
        // negative sign.
        for (double k : signedCurvature(poly))
            signedCurvatures.push_back(sign * k);

        assert(signedCurvatures.size() == vertices.size());
        ScalarField<double> kappa(vertices.size());
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
        ScalarField<double> lengths(vertices.size()); // Actually per edge
        for (size_t i = 0; i < vertices.size(); ++i) {
            auto pt1 = truncateFrom3D<Point2D>(vertices[i]),
                 pt2 = truncateFrom3D<Point2D>(vertices[(i + 1) % vertices.size()]);
            int numAverage = 0;
            double k = 0;
            if (!PeriodicBoundaryMatcher::FaceMembership<2>(pt1, slice.boundingBox()).onAnyFace()) {
                ++numAverage; k += kappa[i];
            }
            if (!PeriodicBoundaryMatcher::FaceMembership<2>(pt2, slice.boundingBox()).onAnyFace()) {
                ++numAverage; k += kappa[(i + 1) % vertices.size()];
            }
            if (numAverage != 0) k /= numAverage;

            // Want to interpolate from upper at k >= c to lower at k = d
            // using function a * 2^(k * b)
            double upper = maxLen / 4.0;
            double lower = minLen;
            double c = -1.0, d = -4;
            // upper = a * 2^(cb)
            // lower = a * 2^(db)
            // upper / lower = 2^((c - d) b) ==> b = log2(u / l) / (c - d)
            // a = upper / 2^(c * b)
            double b = log(upper / lower) / (log(2) * (c - d));
            double a = upper / pow(2, c * b);
            lengths[i] = a * pow(2, b * k);
            lengths[i] = std::min(lengths[i], upper);
            lengths[i] = std::max(lengths[i], lower);
        }
#if DEBUG_OUT
        writer.addField("min_lengths", lengths, DomainType::PER_ELEMENT);
        ScalarField<double> edgeLengths(vertices.size());
        ScalarField<double> isShort(vertices.size());
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
        std::vector<double> &vml = variableMinLens.back();
        for (size_t i = 0; i < lengths.domainSize(); ++i)
            vml[i] = lengths[i];
    }
}
