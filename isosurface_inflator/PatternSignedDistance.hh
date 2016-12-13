////////////////////////////////////////////////////////////////////////////////
// PatternSignedDistance.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computes the signed distance to a pattern represented by the inflation
//      of a WireMesh.
//
//      Blending parameters control the smooth minimum operation's parameter, k,
//      in a nonlinear way:
//          k = Log[2]/s
//      This function is chosen such that the maximum normal shape velocity is
//      unit to match the shape velocity magnitude of the thickness and
//      positional parameters.
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
#include "Joint.hh"
#include "TesselateSpheres.hh"

#include <Future.hh>

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
        m_wireMesh.inflationGraph(params, points, m_edges, thicknesses, m_blendingParams);

        // Vector of edge geometry uses same index as edges
        m_edgeGeometry.clear();
        for (const auto &e : m_edges) {
            m_edgeGeometry.emplace_back(
                    points[e.first],      points[e.second],
                    thicknesses[e.first], thicknesses[e.second]);
        }

        m_adjEdges.resize(points.size());
        for (auto &ae : m_adjEdges) ae.clear();
        for (size_t ei = 0; ei < m_edges.size(); ++ei) {
            const auto &e = m_edges[ei];
            m_adjEdges.at(e.first ).push_back(ei);
            m_adjEdges.at(e.second).push_back(ei);
        }

        // Construct joints at each vertex in the base unit.
        for (size_t u = 0; u < points.size(); ++u) {
            if (m_adjEdges[u].size() < 2) {
                if (WMesh::PatternSymmetry::inBaseUnit(stripAutoDiff(points[u])))
                    throw std::runtime_error("Dangling edge inside base unit");
                // Vertices outside the base unit cell with only one edge
                // incident are not really joints.
                m_jointForVertex.emplace_back(std::unique_ptr<Joint<Real>>());
                continue;
            }
            std::vector<Point3<Real>> centers(1, points[u]);
            std::vector<Real>         radii(1, thicknesses[u]);
            for (size_t ei : m_adjEdges[u]) {
                size_t v_other = m_edges[ei].first == u ? m_edges[ei].second
                                                        : m_edges[ei].first;
                centers.push_back(points[v_other]);
                radii.push_back(thicknesses[v_other]);
            }
            m_jointForVertex.emplace_back(
                Future::make_unique<Joint<Real>>(centers, radii, m_blendingParams[u]));
        }

#if 0
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
                size_t v1 = m_edges[e1].first;
                if (v1 == u) {
                    uIsP1OfE1 = true;
                    v1 = m_edges[e1].second;
                }
                assert(v1 != u);
                Real angleDeficit1 = uIsP1OfE1 ? m_edgeGeometry.at(e1).angleAtP1()
                                               : m_edgeGeometry.at(e1).angleAtP2();
                for (size_t e2 : m_adjEdges[u]) {
                    if (e2 <= e1) continue;
                    bool uIsP1OfE2 = false;
                    size_t v2 = m_edges[e2].first;
                    if (v2 == u) {
                        uIsP1OfE2 = true;
                        v2 = m_edges[e2].second;
                    }
                    assert(v2 != u);

                    Point3<Real> l1 = points[v1] - points[u];
                    Point3<Real> l2 = points[v2] - points[u];
                    l1 /= sqrt(l1.squaredNorm()), l2 /= sqrt(l2.squaredNorm());
                    // get unsigned angle between edges (in [0, Pi])
                    // Real edgeAngle = acos(l1.dot(l2));
                    // std::cout << "ea: " << edgeAngle << std::endl;
                    Real cosTheta = l1.dot(l2);
                    Real sinTheta = sqrt(l1.cross(l2).squaredNorm());
                    // Adept doesn't support atan2...
                    // Subtlety: prevent singularity (nan) in derivative
                    // computation. Inverse trig functions acos and asin aren't
                    // differentiable at +/- 1.0.
                    // We want to use asin when the angle is [0, pi/4), 
                    // acos when the angle is in [pi/4, 3pi/4) and pi + asin
                    // when the angle is in [3pi/4, pi].
                    //
                    // Finally, sinTheta itself is nondifferentiable near 0, so
                    // we explicitly assign 0 to it to avoid nans in automatic
                    // differentiation (will get zero derivative effect)
                    // When sinTheta = 0, we interpret the angle as always
                    // increasing
                    if (sinTheta < 1e-5) sinTheta = 0.0;

                    // First, determine if we're in the left or right quadrant
                    Real edgeAngle;
                    if (cosTheta >= 0.0) {
                        // Right quadrant
                        // Now determine if angle >= pi/4 or < pi/4
                        if (sinTheta >= cosTheta) { edgeAngle = acos(cosTheta); }
                        else                      { edgeAngle = asin(sinTheta); }
                    }
                    else {
                        // Left quadrant
                        // Now determine if angle < 3 pi/4 or >= 3 pi/4
                        if (sinTheta > -cosTheta) { edgeAngle = acos(cosTheta); }
                        else                      { edgeAngle = M_PI + asin(sinTheta); }
                    }

                    // Real initEdgeAngle = edgeAngle;

                    edgeAngle += angleDeficit1;
                    edgeAngle += uIsP1OfE2 ? m_edgeGeometry.at(e2).angleAtP1() : m_edgeGeometry.at(e2).angleAtP2();
                    theta = std::min(theta, edgeAngle);
                    // std::cout << "Vertex " << u
                    //           << " edge pair " << e1 << "," << e2
                    //           << ": initEdgeAngle " << initEdgeAngle
                    //           << ", edgeAngle " << edgeAngle
                    //           << ", angleDeficit1 " << angleDeficit1
                    //           << std::endl;
                }
            }
            Real sinSq4Theta = sin(4.0 * theta);
            // std::cout << "Vertex " << u << " theta: " << theta << std::endl;
            sinSq4Theta *= sinSq4Theta;
            if (theta >= M_PI)               m_vertexSmoothness.push_back(0.0);
            else if (theta > 7 * M_PI / 8.0) m_vertexSmoothness.push_back(sinSq4Theta);
            else                             m_vertexSmoothness.push_back(1.0);
        }

        // for (size_t u = 0; u < points.size(); ++u) {
        //     std::cout << "vertex " << u << " (valence " << m_adjEdges[u].size()
        //          << ") smoothness: " << m_vertexSmoothness[u] << std::endl;
        // }
#endif
        writeDebugSphereMesh("sphere_mesh.msh");
        debugJointAtVertex(4);
    }

    // Additional Real type to support automatic differentiation wrt. p only
    template<typename Real2>
    Real2 signedDistance(Point3<Real2> p) const {
        p = WMesh::PatternSymmetry::mapToBaseUnit(p);
        std::vector<Real2> edgeDists;
        Real2 closestEdgeDist = 1e5;
        size_t closestEdge = m_edgeGeometry.size();
        edgeDists.reserve(m_edgeGeometry.size());
        for (size_t i = 0; i < m_edgeGeometry.size(); ++i) {
            edgeDists.push_back(m_edgeGeometry[i].signedDistance(p));
            if (edgeDists.back() < closestEdgeDist) {
                closestEdgeDist = edgeDists.back();
                closestEdge = i;
            }
        }
        assert(closestEdge < m_edgeGeometry.size());

        std::vector<Real2> jointEdgeDists;
        auto distToVtxJoint = [&](size_t vtx) -> Real2 {
            const auto &joint = m_jointForVertex[vtx];
            if (!joint) {
                // Joints are not created for valence 1 vertices.
                assert(m_adjEdges[vtx].size() == 1);
                return edgeDists[m_adjEdges[vtx][0]];
            }
            jointEdgeDists.clear(), jointEdgeDists.reserve(m_adjEdges[vtx].size());
            for (size_t ei : m_adjEdges[vtx])
                jointEdgeDists.push_back(edgeDists[ei]);
            Real2 s = joint->smoothingAmt(p);
            return SD::exp_smin_reparam_accurate<Real2>(jointEdgeDists, s);
        };

        // Real2 dist = std::min(distToVtxJoint(m_edges[closestEdge].first),
        //                       distToVtxJoint(m_edges[closestEdge].second));

        // Create smoothed union geometry around each vertex and then union
        // together
        Real2 dist = 1e5;
        for (size_t u = 0; u < numVertices(); ++u)
            dist = std::min(dist, distToVtxJoint(u));
#if 0
        for (size_t u = 0; u < numVertices(); ++u) {
            // Vertex smoothness is in [0, 1],
            //      1.0: full smoothness (m_blendingParams(u))
            //      0.0: "no" smoothness (1/256.0)
            // "smoothness" is computed in multiple steps to work around a bug
            // in Eigen AutoDiff's make_coherent for expression templates.
            Real2 smoothness = 1/256.0 + (1.0 - 1/256.0) * vertexSmoothness(u);
            smoothness *= m_blendingParams.at(u);
            // Transition to precise union for extremely low smoothing
            if (smoothness > 1/256.0) {
                // inlined exp_smin_reparam
                Real2 k = 1 / smoothness;
                Real2 smin = 0;
                for (size_t ei : m_adjEdges[u])
                    smin += exp(-k * edgeDists[ei]);
                dist = std::min<Real2>(dist, -log(smin) * smoothness);
            }
            else {
                for (size_t ei : m_adjEdges[u])
                    dist = std::min<Real2>(dist, edgeDists.at(ei));
            }
        }
#endif

        assert(!std::isnan(stripAutoDiff(dist)));
        assert(!std::isinf(stripAutoDiff(dist)));

        return dist;
    }

    // Accelerated version of "signedDistance(p) <= 0"
    bool isInside(Point3<Real> p) const {
        p = WMesh::PatternSymmetry::mapToBaseUnit(p);
        std::vector<Real> edgeDists;
        edgeDists.reserve(m_edgeGeometry.size());
        // Definitely inside if we're inside one of the edges: assumes blending
        // is additive.
        // Note: possibly could be sped up by implementing cheap isInside for
        // edge geometry and bailing early if inside?
        Real closestEdgeDist = 1e5;
        size_t closestEdge = m_edgeGeometry.size();
        for (size_t i = 0; i < m_edgeGeometry.size(); ++i) {
            edgeDists.push_back(m_edgeGeometry[i].signedDistance(p));
            if (edgeDists.back() < closestEdgeDist) {
                closestEdgeDist = edgeDists.back();
                closestEdge = i;
            }
        }
        assert(closestEdge < m_edgeGeometry.size());

        std::vector<Real> jointEdgeDists;
        auto distToVtxJoint = [&](size_t vtx) {
            const auto &joint = m_jointForVertex[vtx];
            if (!joint) {
                // Joints are not created for valence 1 vertices.
                assert(m_adjEdges[vtx].size() == 1);
                return edgeDists[m_adjEdges[vtx][0]];
            }
            jointEdgeDists.clear(), jointEdgeDists.reserve(m_adjEdges[vtx].size());
            for (size_t ei : m_adjEdges[vtx])
                jointEdgeDists.push_back(edgeDists[ei]);
            Real s = joint->smoothingAmt(p);
            return SD::exp_smin_reparam_accurate<Real>(jointEdgeDists, s);
        };

        return std::min(distToVtxJoint(m_edges[closestEdge].first),
                        distToVtxJoint(m_edges[closestEdge].second)) <= 0;

        // // see signedDistance(p) above
        // for (size_t u = 0; u < numVertices(); ++u)
        //     if (distToVtxJoint(u) <= 0) return true;
#if 0
        for (size_t u = 0; u < numVertices(); ++u) {
            Real smoothness = m_blendingParams[u] * (1/256.0 + (1.0 - 1/256.0) * m_vertexSmoothness[u]);
            if (smoothness > 1/256.0) {
                const Real k = 1.0 / smoothness;
                Real res = 0;
                for (size_t ei : m_adjEdges[u])
                    res += exp(-k * edgeDists[ei]);
                if (res > 1) return true;
            }
            else { continue; } // Note: sharp joints only contain p if one of the edges does (can't happen)
        }
#endif

        return false;
    }

    // Debug smoothing modulation field/smoothing amount
    template<typename Real2>
    Real2 smoothnessAtEvalPt(Point3<Real2> p) const {
        p = WMesh::PatternSymmetry::mapToBaseUnit(p);
        std::vector<Real2> edgeDists;
        Real2 closestEdgeDist = 1e5;
        size_t closestEdge = m_edgeGeometry.size();
        edgeDists.reserve(m_edgeGeometry.size());
        for (size_t i = 0; i < m_edgeGeometry.size(); ++i) {
            edgeDists.push_back(m_edgeGeometry[i].signedDistance(p));
            if (edgeDists.back() < closestEdgeDist) {
                closestEdgeDist = edgeDists.back();
                closestEdge = i;
            }
        }
        assert(closestEdge < m_edgeGeometry.size());

        // Create smoothed union geometry around each vertex and then union
        // together
        Real2 dist = 1e5;
        Real2 smoothness = -1.0;
        std::vector<Real2> jointEdgeDists;
        for (size_t u = 0; u < numVertices(); ++u) {
            const auto &joint = m_jointForVertex[u];
            Real2 jdist;
            Real2 s;
            if (!joint) {
                // Joints are not created for valence 1 vertices.
                assert(m_adjEdges[u].size() == 1);
                jdist = edgeDists[m_adjEdges[u][0]]; 
                s = m_blendingParams[u];
            }
            else {
                jointEdgeDists.clear(), jointEdgeDists.reserve(m_adjEdges[u].size());
                for (size_t ei : m_adjEdges[u])
                    jointEdgeDists.push_back(edgeDists[ei]);
                s = joint->smoothingAmt(p);
                jdist = SD::exp_smin_reparam_accurate<Real2>(jointEdgeDists, s);
            }
            if (jdist < dist) {
                dist = jdist;
                smoothness = s;
            }
        }

        return smoothness;
    }

    // Representative cell bounding box (region to be meshed)
    BBox<Point3<Real>> boundingBox() const { return m_bbox; }
    // Choose a different box region to be meshed instead of the default
    // representative cell region (for debugging purposes)
    void setBoundingBox(const BBox<Point3<Real>> &bb) { m_bbox = bb; }

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
            for (size_t d = 0; d < 3; ++d) {
                boxCorner[d] = (corner & (1 << d)) ? bbox.maxCorner[d] : bbox.minCorner[d];
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
    size_t numVertices() const { return m_jointForVertex.size(); }

    // Will differ from numVertices() since valence-1 vertices (outside the
    // meshing cell) do not have joints.
    size_t numJoints() const {
        size_t count = 0;
        for (auto &joint : m_jointForVertex)
            if (joint) ++count;
        return count;
    }

    // Write a mesh consisting of a sphere at each joint tagged with the joint's
    // vertex index and its blending parameters.
    void writeDebugSphereMesh(const std::string &path) {
        const size_t numSpherePoints = 1000;

        std::vector<Point3<double>> centers;
        std::vector<double>         radii;
        std::vector<size_t>         vertexForJoint;
        for (size_t i = 0; i < m_jointForVertex.size(); ++i) {
            auto &joint = m_jointForVertex[i];
            if (joint) {
                centers.push_back(stripAutoDiff(joint->c1()));
                radii  .push_back(stripAutoDiff(joint->r1()));
                vertexForJoint.push_back(i);
            }
        }

        std::vector<MeshIO::IOVertex > outVertices;
        std::vector<MeshIO::IOElement> outElements;
        std::vector<size_t> sphereIdx;
        tesselateSpheres(numSpherePoints, centers, radii, outVertices, outElements, sphereIdx);

        ScalarField<double> vertexIndex(outVertices.size()),
                            blendingParam(outVertices.size());
        for (size_t i = 0; i < outVertices.size(); ++i) {
            vertexIndex[i]   = vertexForJoint.at(sphereIdx[i]);
            blendingParam[i] = stripAutoDiff(m_blendingParams[vertexForJoint.at(sphereIdx[i])]);
        }

        MSHFieldWriter writer(path, outVertices, outElements);
        writer.addField("vtx_index",     vertexIndex, DomainType::PER_NODE);
        writer.addField("blend_param", blendingParam, DomainType::PER_NODE);
    }

    void debugJointAtVertex(size_t vertexIndex) const {
        const size_t numSpherePoints = 1000;
        auto &joint = m_jointForVertex[vertexIndex];
        if (!joint) throw std::runtime_error("No joint at vertex " + std::to_string(vertexIndex));

        const auto &hull = joint->blendingHull();
        std::vector<Point3<double>> centers;
        std::vector<double>         radii;
        for (auto &pt : hull.sphereCenters()) centers.push_back(stripAutoDiff(pt));
        for (auto  &r : hull.sphereRadii  ()) radii.push_back(stripAutoDiff(r));

        std::vector<MeshIO::IOVertex > outVertices;
        std::vector<MeshIO::IOElement> outElements;
        std::vector<size_t> sphereIdx;
        tesselateSpheres(numSpherePoints, centers, radii, outVertices, outElements, sphereIdx);
        std::string name = "joint" + std::to_string(vertexIndex);
        MSHFieldWriter writer(name + ".msh", outVertices, outElements);

        std::ofstream outFile(name + ".txt");
        outFile << std::setprecision(19);
        outFile << "Sphere centers:" << std::endl;
        for (auto &pt : centers) outFile << "    {" << pt.transpose() << "}" << std::endl;
        outFile << "Sphere radii:" << "{";
        for (double r : radii) outFile << r << ", ";
        outFile << "}" << std::endl;
    }

private:
    const WMesh &m_wireMesh;

    // Bounding box for the meshing cell. Defaults to the representative cell
    // for the symmetry type, but can be changed manually for debugging
    // purposes.
    BBox<Point3<Real>> m_bbox = WMesh::PatternSymmetry::template representativeMeshCell<Real>();

    std::vector<typename WMesh::Edge> m_edges; // vertex index pairs for each edge
    std::vector<SD::Primitives::InflatedEdge<Real>> m_edgeGeometry;
    // Quantity between 0 and 1 saying how much to smooth intersections between
    // edges incident on a particular vertex
    // Since a given edge pair only has (at most) one vertex in common,
    // specifying this per-vertex makes sense.
    std::vector<Real>                m_vertexSmoothness;
    std::vector<Real>                m_blendingParams;
    // Edges adjacent to a particular vertex.
    std::vector<std::vector<size_t>> m_adjEdges;

    // Joint vertex indices and their geometry
    // (Not every vertex is a joint; there are dangling edges extending outside
    //  the base unit.)
    std::vector<std::unique_ptr<Joint<Real>>> m_jointForVertex;
};

#endif /* end of include guard: PATTERNSIGNEDDISTANCE_HH */
