////////////////////////////////////////////////////////////////////////////////
// PatternSignedDistance.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computes the signed distance to a pattern represented by the inflation
//      of a WireMesh.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/23/2015 14:58:31
////////////////////////////////////////////////////////////////////////////////
#ifndef PATTERNSIGNEDDISTANCE_HH
#define PATTERNSIGNEDDISTANCE_HH
#include "WireMesh.hh"
#include "AutomaticDifferentiation.hh"
#include "SignedDistance.hh"

template<typename _Real, class WMesh>
class PatternSignedDistance {
public:
    typedef _Real Real;
    PatternSignedDistance(const WMesh &wireMesh) : m_wireMesh(wireMesh) {
        static_assert(WMesh::thicknessType == ThicknessType::Vertex,
                "Only per-vertex thicknesses are currently supported");
    }

    size_t numParams() const { return m_wireMesh.numParams(); }
    size_t numThicknessParams() const { return m_wireMesh.numThicknessParams(); }
    size_t numPositionParams() const { return m_wireMesh.numPositionParams(); }
    size_t numBlendingParams() const { return m_wireMesh.numBlendingParams(); }
    void setParameters(const std::vector<Real> &params) {
        std::vector<Real> thicknesses;
        std::vector<Point3<Real>> points;
        std::vector<typename WMesh::Edge> edges;
        m_wireMesh.inflationGraph(params, points, edges, thicknesses, m_blendingParams);

        // Vector of edge geometry uses same index as edges
        m_edgeGeometry.clear();
        for (const auto &e : edges) {
            m_edgeGeometry.emplace_back(
                    points[e.first],      points[e.second],
                    thicknesses[e.first], thicknesses[e.second]);
        }

        m_adjEdges.resize(points.size());
        for (auto &ae : m_adjEdges) ae.clear();
        for (size_t ei = 0; ei < edges.size(); ++ei) {
            const auto &e = edges[ei];
            m_adjEdges.at(e.first ).push_back(ei);
            m_adjEdges.at(e.second).push_back(ei);
        }

        // Compute vertex smoothness:
        // Vertices with intersecting edges are smoothed. This smoothing is
        // ramped up from 0.0 to 1.0 as the minimum incident angle, "theta"
        // shrinks from Pi to Pi/2:
        // smoothness = 0             if theta >= Pi
        //              sin^2(theta)  if Pi/2 < theta < Pi
        //              1.0           if theta < Pi/2
        m_vertexSmoothness.clear(), m_vertexSmoothness.reserve(points.size());
        for (size_t u = 0; u < points.size(); ++u) {
            // Min angle over all pairs of edges
            Real theta = 2 * M_PI;
            for (size_t e1 : m_adjEdges[u]) {
                // Get other vertex of e1, and determine if it is that edge's
                // "p1" or "p2" endpoint.
                bool uIsP1OfE1 = false;
                size_t v1 = edges[e1].first;
                if (v1 == u) {
                    uIsP1OfE1 = true;
                    v1 = edges[e1].second;
                }
                assert(v1 != u);
                Real angleDeficit1 = uIsP1OfE1 ? m_edgeGeometry.at(e1).angleAtP1()
                                               : m_edgeGeometry.at(e1).angleAtP2();
                for (size_t e2 : m_adjEdges[u]) {
                    if (e2 <= e1) continue;
                    bool uIsP1OfE2 = false;
                    size_t v2 = edges[e2].first;
                    if (v2 == u) {
                        uIsP1OfE2 = true;
                        v2 = edges[e2].second;
                    }
                    assert(v2 != u);

                    Point3<Real> l1 = points[v1] - points[u];
                    Point3<Real> l2 = points[v2] - points[u];
                    l1 /= sqrt(l1.squaredNorm()), l2 /= sqrt(l2.squaredNorm());
                    // get angle between edges (in [0, Pi])
                    Real edgeAngle = acos(l1.dot(l2));
                    edgeAngle -= angleDeficit1;
                    edgeAngle -= uIsP1OfE2 ? m_edgeGeometry.at(e2).angleAtP1() : m_edgeGeometry.at(e2).angleAtP2();
                    theta = std::min(theta, edgeAngle);
                }
            }
            Real sinSqTheta = sin(theta);
            // cout << "Vertex " << u << " theta: " << theta << endl;
            sinSqTheta *= sinSqTheta;
            if (theta >= M_PI)           m_vertexSmoothness.push_back(0.0);
            else if (theta > M_PI / 2.0) m_vertexSmoothness.push_back(sinSqTheta);
            else                         m_vertexSmoothness.push_back(1.0);
        }

        // for (size_t u = 0; u < points.size(); ++u) {
        //     cout << "vertex " << u << " (valence " << m_adjEdges[u].size()
        //          << ") smoothness: " << m_vertexSmoothness[u] << endl;
        // }
    }

    // Additional Real type to support automatic differentiation wrt. p only
    template<typename Real2>
    Real2 signedDistance(Point3<Real2> p) const {
        p = WMesh::PatternSymmetry::mapToBaseUnit(p);
        std::vector<Real2> edgeDists;
        edgeDists.reserve(m_edgeGeometry.size());
        for (const auto &c : m_edgeGeometry)
            edgeDists.push_back(c.signedDistance(p));

        // for (Real2 ed : edgeDists) {
        //     std::cout << "Edge distance: " << ed << std::endl;
        // }

        // Create smoothed union geometry around each vertex and then union
        // together
        std::vector<Real2> unionedDist;
        std::vector<Real2> incidentDists;
        Real2 dist = 1e5;
        for (size_t u = 0; u < numVertices(); ++u) {
            incidentDists.clear();
            for (size_t ei : m_adjEdges[u])
                incidentDists.push_back(edgeDists.at(ei));
            dist = std::min(dist, SD::mix(SD::min(incidentDists), SD::exp_smin(incidentDists, Real2(m_blendingParams.at(u))),
                                 Real2(vertexSmoothness(u))));
        }

        return dist;
    }

    bool isInside(Point3<Real> &p) const {
        return signedDistance(p) <= 0;
    }

    // Representative cell bounding box.
    static BBox<Point3<Real>> boundingBox() {
        return WMesh::PatternSymmetry::template representativeMeshCell<Real>();
    }

    // Sphere bounding the representative mesh cell, needed for CGAL meshing.
    // Note: CGAL requires the bounding sphere center to lie inside the object.
    template<typename Real>
    void boundingSphere(Point3<Real> &c, Real &r) const {
        auto bbox = boundingBox();
        auto boxCenter = bbox.center();
        // CGAL requires our bounding sphere center to lie within the object, so we
        // find the closest point on the medial axis to the bounding box center.
        c = closestMedialAxisPoint(boxCenter);
        // Determine a radius large enough for the bounding box to fit within.
        Point3<Real> boxCorner;
        r = 0.0;
        for (size_t corner = 0; corner < 8; ++corner) {
            for (size_t c = 0; c < 3; ++c) {
                boxCorner[c] = corner & (1 << c) ? bbox.maxCorner[c] : bbox.minCorner[c];
            }
            r = std::max(r, (c - boxCorner).norm());
        }
        r += 1e-2; // Add some padding to the bounding sphere.
    }

    Point3<Real> closestMedialAxisPoint(const Point3<Real> &p) const {
        Real dist = std::numeric_limits<Real>::max();
        Point3<Real> closestPoint(Point3<Real>::Zero());
        for (const auto &edge : m_edgeGeometry) {
            auto c = edge.closestMedialAxisPoint(p);
            Real candidateDist = (c - p).squaredNorm();
            if (candidateDist < dist) {
                dist = candidateDist;
                closestPoint = c;
            }
        }
        return closestPoint;
    }

    Real vertexSmoothness(size_t v) const { return m_vertexSmoothness.at(v); }
    size_t numVertices() const { return m_vertexSmoothness.size(); }

private:
    const WMesh &m_wireMesh;

    std::vector<SD::Primitives::InflatedEdge<Real>> m_edgeGeometry;
    // Quantity between 0 and 1 saying how much to smooth intersections between
    // edges incident on a particular vertex
    // Since a given edge pair only has (at most) one vertex in common,
    // specifying this per-vertex makes sense.
    std::vector<Real>                m_vertexSmoothness;
    std::vector<Real>                m_blendingParams;
    // Edges adjacent to a particular vertex.
    std::vector<std::vector<size_t>> m_adjEdges;
};

#endif /* end of include guard: PATTERNSIGNEDDISTANCE_HH */
