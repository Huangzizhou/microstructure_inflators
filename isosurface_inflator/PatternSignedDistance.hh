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

#include "SignedDistanceRegion.hh"

#include <Future.hh>

#define VERTEX_SMOOTHNESS_MODULATION 1
#define DISCONTINUITY_AVOIDING_CREASE_AVOIDANCE 1

template<typename _Real, class WMesh>
class PatternSignedDistance : public SignedDistanceRegion<3> {
public:
    using Real = _Real;
    PatternSignedDistance(const WMesh &wireMesh) : m_wireMesh(wireMesh) {
        if (wireMesh.thicknessType() != ThicknessType::Vertex)
            throw std::runtime_error("Only Per-Vertex Thickness Supported");
    }

    // Always support double type for compatibility with SignedDistanceRegion
    virtual double signedDistance(const Point3D &p) const override { return stripAutoDiff(m_signedDistanceImpl(autodiffCast<Point3<Real>>(p))); }

    // Also support automatic differentiation types
    template<typename Real2, bool DebugDerivatives = false>
    Real2 signedDistance(const Point3<Real2> &p) const { return m_signedDistanceImpl(p); }

    size_t numParams() const { return m_wireMesh.numParams(); }
    size_t numThicknessParams() const { return m_wireMesh.numThicknessParams(); }
    size_t numPositionParams() const { return m_wireMesh.numPositionParams(); }
    size_t numBlendingParams() const { return m_wireMesh.numBlendingParams(); }

    // Accelerated version of "signedDistance(p) <= 0"
    virtual bool isInside(const Point3D &p) const override {
        return m_isInsideImpl(p);
    }

    void setParameters(const std::vector<Real> &params,
                       JointBlendMode blendMode = JointBlendMode::HULL) {
        // Clear all existing state.
        m_edgeGeometry.clear(), m_jointForVertex.clear(),
        m_vertexSmoothness.clear(), m_adjEdges.clear();

        std::vector<Real> thicknesses;
        std::vector<Point3<Real>> points;
        m_wireMesh.inflationGraph(params, points, m_edges, thicknesses, m_blendingParams);

        // Vector of edge geometry uses same index as edges
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
                if (WMesh::PatternSymmetry::inBaseUnit(stripAutoDiff(points[u]))) {
                    std::cerr << "WARNING: dangling edge inside the base unit at ["
                              << points[u].transpose()
                              << "] (neighboring cell's edge is probably protruding inside base cell)"
                              << std::endl;
                }
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
                Future::make_unique<Joint<Real>>(centers, radii, m_blendingParams[u], blendMode));
        }

#if VERTEX_SMOOTHNESS_MODULATION
        // Compute vertex smoothness:
        // Vertices with intersecting edges are smoothed. This smoothing is
        // ramped up from 0.0 to 1.0 as the minimum incident angle, "theta"
        // shrinks from Pi to Pi/2:
        // smoothness = 0             if theta >= Pi
        //              sin^2(theta)  if Pi/2 < theta < Pi
        //              1.0           if theta < Pi/2
        m_vertexSmoothness.reserve(points.size());
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
#if 0
        writeDebugSphereMesh("sphere_mesh.msh");
        debugJointAtVertex(4);
#endif
#if 0
        // For debugging disconnected component in 0746:
        if (!isAutodiffType(Real(0))) {
            signedDistance<Real, true>(Point3<Real>(0.28516, 0.01, 0.998));
            signedDistance<Real, true>(Point3<Real>(0.28516, 0.005, 0.998));
            signedDistance<Real, true>(Point3<Real>(0.28516, 0.0025, 0.998));
            signedDistance<Real, true>(Point3<Real>(0.28516, 0, 0.998));
        }
#endif
    }

    // Distance to both the smoothed and hard-unioned versions of a joint.
    template<typename Real2>
    struct JointDists {
        JointDists(Real2 sd, Real2 hd) : smooth(sd), hard(hd) { }
        JointDists() { } // Uninitialized!
        Real2 smooth, hard;
        static JointDists largest() {
            return JointDists{safe_numeric_limits<Real2>::max(), safe_numeric_limits<Real2>::max()};
        }
        // The geometry is determined by the smooth distance, so this should be
        // used to sort/compare distance pairs
        bool operator<(const JointDists<Real2> &b) const { return smooth < b.smooth; }
    };

    // Determine distance from "p" to the joint at vertex "vtx" whose edges are
    // each at signed distance "edgeDists".
    // "jointEdgeDists" is used as scratch space to prevent allocation
    // Returns
    //      distance to smoothed joint (first)
    //      distance to hard-unioned joint (second)
    template<typename Real2, bool DebugOutput = false>
    JointDists<Real2>
    distToVtxJoint(size_t vtx, const Point3<Real2> &p,
                   const std::vector<Real2> &edgeDists,
                   std::vector<Real2> &jointEdgeDists) const {
        const auto &joint = m_jointForVertex[vtx];

        if (!joint) {
            // Joints are not created for valence 1 vertices.
            assert(m_adjEdges[vtx].size() == 1);
            return JointDists<Real2>(edgeDists[m_adjEdges[vtx][0]],
                                     edgeDists[m_adjEdges[vtx][0]]);
        }
        jointEdgeDists.clear(), jointEdgeDists.reserve(m_adjEdges[vtx].size());
        Real2 hardUnionedDist = safe_numeric_limits<Real2>::max();
        for (size_t ei : m_adjEdges[vtx]) {
            jointEdgeDists.push_back(edgeDists[ei]);
            hardUnionedDist = std::min<Real2>(hardUnionedDist, edgeDists[ei]);
        }
        if (DebugOutput) {
            std::cerr << "hardUnionedDist derivatives:";
            reportDerivatives(std::cerr, hardUnionedDist);
            std::cerr << std::endl;
        }
        Real2 s = joint->template smoothingAmt<Real2, DebugOutput>(p);
        if (VERTEX_SMOOTHNESS_MODULATION)
            s *= m_vertexSmoothness[vtx];
        if (DebugOutput) {
            std::cerr << "js derivatives:"; reportDerivatives(std::cerr,  joint->smoothingAmt(p)); std::cerr << std::endl;
            std::cerr << "vs derivatives:"; reportDerivatives(std::cerr, m_vertexSmoothness[vtx]); std::cerr << std::endl;
        }
        return JointDists<Real2>(SD::exp_smin_reparam_accurate<Real2>(jointEdgeDists, s),
                                 hardUnionedDist);
    }

    // InsideOutsideAccelerate: do we just need an inside/outside query? If so,
    // we can implement some optimizations
    template<typename Real2, bool InsideOutsideAccelerate = false, bool DebugDerivatives = false>
    Real2
    combinedJointDistances(const Point3<Real2> &p,
                           const std::vector<Real2> &edgeDists, size_t /* closestEdge */) const {
        std::vector<Real2> jointEdgeDists;
        const double maxOverlapSmoothingAmt = 0.02;

        // {
        //     Real2 dist = 1e5;
        //     for (size_t vtx = 0, i = 0; vtx < numVertices(); ++vtx)
        //         dist = std::min(dist, distToVtxJoint(vtx, p, edgeDists, jointEdgeDists).smooth);
        //     return dist;
        // }
#if 1
        // Note: this computation is made slow by needing to compute signed
        // distances to the blending region for each joint considered. To
        // accelerate things, we reduce the number of joints we smooth.
        // We just need to make sure the two closest joints to p, in terms of
        // smoothed signed distance, are considered. It is difficult to
        // predict which two joints these are from edge distances alone as
        // the smoothing can change which joint is closest.
        //
        // For now, we make the assumption that only the closest N of the joints
        // in terms of hard-unioned distance are candidates for the closest two
        // smoothed joints.

        // Compute hard-unioned distance to each joint and determine Nth closest
        std::vector<double> hard_distance(numVertices(), safe_numeric_limits<double>::max());
        for (size_t vtx = 0; vtx < numVertices(); ++vtx) {
            for (size_t ei : m_adjEdges[vtx])
                hard_distance[vtx] = std::min<double>(hard_distance[vtx],
                                                      stripAutoDiff(edgeDists[ei]));
        }

        double candidateDistThreshold;
        const size_t MAX_CANDIDATES = 300; // conservative upper bound for array allocation.
        size_t numCandidates = 0;
        {
            size_t requestedCandidates = std::min<size_t>(20, hard_distance.size() - 1);
            std::vector<double> hd_copy(hard_distance);
            std::nth_element(hd_copy.begin(), hd_copy.begin() + requestedCandidates,
                             hd_copy.end());
            candidateDistThreshold = hd_copy[requestedCandidates];
            for (double hd : hard_distance)
                if (hd <= candidateDistThreshold) ++numCandidates;
        }
        if (numCandidates > MAX_CANDIDATES) {
            std::cerr << "numCandidates: " << numCandidates << std::endl;
            assert(numCandidates <= MAX_CANDIDATES);
        }

        static_assert(!(InsideOutsideAccelerate && isAutodiffType<Real2>()),
                      "The inside-outside test is non-differentiable");
        if (InsideOutsideAccelerate) {
            double conservativeSMin = 0.0;
            for (size_t vtx = 0; vtx < numVertices(); ++vtx) {
                if (hard_distance[vtx] > candidateDistThreshold) continue; // prune far joints
                double joint_smin = 0;
                // m_vertexSmoothness[vtx] could scale down smoothing--but for
                // robustness (to avoid small smoothing params) we
                // conservatively assume m_vertexSmoothness == 1.
                double s = stripAutoDiff(m_blendingParams[vtx]);
                double k = 1.0 / s;
                for (size_t ei : m_adjEdges[vtx])
                    joint_smin += exp(-k * stripAutoDiff(edgeDists[ei]));
                // joint_smin = -s * log(joint_smin);
                //
                // Individual joint_smin values are then smin-ed together with
                // smoothing maxOverlapSmoothingAmt:
                // exp(-(-s * log(joint_smin)) / maxOverlapSmoothingAmt)
                // = exp(log(joint_smin))^(s/maxOverlapSmoothingAmt)
                // = joint_smin^(s/maxOverlapSmoothingAmt)
                conservativeSMin += pow(joint_smin,
                                        s / maxOverlapSmoothingAmt);
                // conservativeSMin += exp(-joint_smin / maxOverlapSmoothingAmt);
            }
            // conservatively inside if
            // -maxOverlapSmoothingAmt * log(conservativeSMin) <= 0
            // <==> log(conservativeSMin) >= 0
            // <==> conservativeSMin >= 1
            //  ==> we are definitely outside if conservativeSMin < 1
            if (conservativeSMin < 1)
                return 1.0;
        }

        // Compute both smoothed and hard-unioned distances to the two closest
        // smoothed joints. These are the joints that could possibly overlap to
        // form a hard crease.
        std::array<JointDists<Real2>, MAX_CANDIDATES> candidateJDists;
        JointDists<Real2> closestJDist, secondClosestJDist;
        closestJDist = secondClosestJDist = JointDists<Real2>::largest();
        size_t  c_idx = safe_numeric_limits<size_t>::max(),
               sc_idx = safe_numeric_limits<size_t>::max();
        for (size_t vtx = 0, i = 0; vtx < numVertices(); ++vtx) {
            if (hard_distance[vtx] > candidateDistThreshold) continue; // prune out the far joints
            candidateJDists[i] = distToVtxJoint(vtx, p, edgeDists, jointEdgeDists);
            const JointDists<Real2> &d = candidateJDists[i];
            if (d < secondClosestJDist) {
                if (d < closestJDist) {
                    secondClosestJDist = closestJDist;
                    closestJDist = d;
                    sc_idx = c_idx;
                    c_idx = vtx;
                }
                else { secondClosestJDist = d; sc_idx = vtx; }
            }
            ++i;
        }

#if DISCONTINUITY_AVOIDING_CREASE_AVOIDANCE
        // Creases can form in the hard-union of all joint geometries when the
        // blending regions of the two joints incident an edge overlap. We avoid
        // this by detecting such region overlaps and applying a smooth union.

        // A point is in the overlap of the blending regions when the two
        // closest joints each has a smooth-union distance differing from its
        // hard-union distance. We thus we smooth to a degree determined by
        // the geometric mean of these differences.

        // Note: if we only consider the two closest joints in computing this
        // smoothing amount, the amount can change abruptly in degenerate cases
        // where many joints are roughly the same distance from the evaluation
        // point. This will case a bad surface discontinuity.

        // To avoid this problem, we take a weighted geometric mean of all
        // differences instead of just the two closest joints' differences.
        // We choose the weighting based on the distance and have it fall off
        // rapidly after the first two closest.
        std::array<Real2, MAX_CANDIDATES> smoothEffects, weights; // Large enough not to matter
        for (size_t i = 0; i < numCandidates; ++i) {
            // Interpolate from 1.0 to 0
            Real2 dist = candidateJDists[i].smooth - secondClosestJDist.smooth;
            // Real2 posDist = std::max<Real2>(dist, 0);
            // weights[i] = 0.5 * (1.0 - tanh(posDist * 100 - 3.0));
            weights[i] = 0.5 * (1.0 - tanh(dist * 50 - 1.0));
            smoothEffects[i] = candidateJDists[i].hard - candidateJDists[i].smooth;
            assert(smoothEffects[i] >= -1e-9);
            smoothEffects[i] = std::max<Real2>(smoothEffects[i], 0.0);
        }
        Real2 weightedAvgOfGeometricMeanSq = 0;
        Real2 totalWeight = 0;
        for (size_t i = 0; i < numCandidates; ++i) {
            for (size_t j = i + 1; j < numCandidates; ++j) {
                weightedAvgOfGeometricMeanSq += weights[i] * weights[j] * smoothEffects[i] * smoothEffects[j];
                totalWeight += weights[i] * weights[j];
            }
        }
        weightedAvgOfGeometricMeanSq /= totalWeight;
        Real2 overlapSmoothAmt = maxOverlapSmoothingAmt * tanh(1000.0 * weightedAvgOfGeometricMeanSq);
#else
        Real2 smoothEffect1 =       closestJDist.hard -       closestJDist.smooth;
        Real2 smoothEffect2 = secondClosestJDist.hard - secondClosestJDist.smooth;
        // Smoothing is additive
        assert(smoothEffect1 >= -1e-9);
        assert(smoothEffect2 >= -1e-9);
        smoothEffect1 = std::max<Real2>(smoothEffect1, 0.0);
        smoothEffect2 = std::max<Real2>(smoothEffect2, 0.0);
        // Choose smoothing based on the gemoetric mean of the differences
        // between distances to smooth and hard-unioned geometry.
        Real2 meanSmoothEffectSq = smoothEffect1 * smoothEffect2;

        // We want to interpolate smoothly from 0 when meanSmoothEffectSq = 0 to ~0.1
        Real2 overlapSmoothAmt = maxOverlapSmoothingAmt * tanh(1000.0 * meanSmoothEffectSq);
#endif

        Real2 dist = SD::exp_smin_reparam_accurate(closestJDist.smooth,
                                             secondClosestJDist.smooth, overlapSmoothAmt);
        if (std::isnan(stripAutoDiff(dist)) || std::isinf(stripAutoDiff(dist))) {
            std::cerr.precision(19);
            std::cerr << "closestJDist.smooth: "       << closestJDist.smooth       << std::endl;
            std::cerr << "secondClosestJDist.smooth: " << secondClosestJDist.smooth << std::endl;
            std::cerr << "closestJDist.hard: "         << closestJDist.hard         << std::endl;
            std::cerr << "secondClosestJDist.hard: "   << secondClosestJDist.hard   << std::endl;
            std::cerr << "overlapSmoothAmt: "          << overlapSmoothAmt          << std::endl;
            std::cerr << "numCandidates: " << numCandidates << std::endl;

            std::cerr << std::endl;
            std::cerr << "candidateDistThreshold: " << candidateDistThreshold << std::endl;
            std::cerr << "hard_distance: ";
            for (Real2 hd : hard_distance) {
                std::cerr << "\t" << hd;
            }
            std::cerr << std::endl;
        }

        // dist = closestJDist.smooth;
        if (hasInvalidDerivatives(dist) || DebugDerivatives) {
            std::cerr << "                     dist:"; std::cerr <<                      dist << std::endl;
            std::cerr << "      closestJDist.smooth:"; std::cerr <<       closestJDist.smooth << std::endl;
            std::cerr << "secondClosestJDist.smooth:"; std::cerr << secondClosestJDist.smooth << std::endl;
            std::cerr << "         overlapSmoothAmt:"; std::cerr <<          overlapSmoothAmt << std::endl;

            const auto & cJoint = m_jointForVertex[c_idx];
            const auto &scJoint = m_jointForVertex[sc_idx];
            Real2  cSmoothness =  cJoint->smoothingAmt(p);
            Real2 scSmoothness = scJoint->smoothingAmt(p);
            if (VERTEX_SMOOTHNESS_MODULATION) {
                 cSmoothness *= m_vertexSmoothness[ c_idx];
                scSmoothness *= m_vertexSmoothness[sc_idx];
            }
            std::cerr << "      closestJDist smoothing amt:"; std::cerr <<  cSmoothness << std::endl;
            std::cerr << "secondClosestJDist smoothing amt:"; std::cerr << scSmoothness << std::endl;

            if (hasInvalidDerivatives(dist))
                std::cerr << "Invalid derivatives computed in combinedJointDistances evaluation" << std::endl;
            std::cerr << "dist:"                     ; reportDerivatives(std::cerr,                      dist); std::cerr << std::endl; std::cerr << std::endl;
            std::cerr << "closestJDist.smooth:"      ; reportDerivatives(std::cerr,       closestJDist.smooth); std::cerr << std::endl; std::cerr << std::endl;
            std::cerr << "secondClosestJDist.smooth:"; reportDerivatives(std::cerr, secondClosestJDist.smooth); std::cerr << std::endl; std::cerr << std::endl;
            std::cerr << "overlapSmoothAmt:"         ; reportDerivatives(std::cerr,          overlapSmoothAmt); std::cerr << std::endl; std::cerr << std::endl;

            std::cerr << "Hard derivatives:" << std::endl;
            std::cerr << "closestJDist.hard:"        ; reportDerivatives(std::cerr,       closestJDist.hard); std::cerr << std::endl;
            std::cerr << "secondClosestJDist.hard:"  ; reportDerivatives(std::cerr, secondClosestJDist.hard); std::cerr << std::endl;
            std::cerr << std::endl;

            for (size_t ei : m_adjEdges[c_idx]) {
                size_t a = m_edges[ei].first, b = m_edges[ei].second;
                std::cerr << "Edge (" << a << ", " << b << ") dist " << edgeDists[ei] << " derivatives:";
                reportDerivatives(std::cerr, edgeDists[ei]);
                std::cerr << std::endl;
            }
            std::cerr << std::endl;

            distToVtxJoint<Real2, true>(c_idx, p, edgeDists, jointEdgeDists);

            if (hasInvalidDerivatives(dist))
                throw std::runtime_error("Invalid derivatives computed in combinedJointDistances evaluation");
        }
#else
        Real2 dist = std::min(distToVtxJoint(m_edges[closestEdge].first, p, edgeDists, jointEdgeDists).smooth,
                              distToVtxJoint(m_edges[closestEdge].second, p, edgeDists, jointEdgeDists).smooth);
#endif
        return dist;
    }

    // Debug smoothing modulation field/smoothing amount
    template<typename Real2>
    std::pair<Real2, size_t> smoothnessAndClosestVtx(Point3<Real2> p) const {
        p = WMesh::PatternSymmetry::mapToBaseUnit(p);
        std::vector<Real2> edgeDists;
        edgeDists.reserve(m_edgeGeometry.size());
        for (size_t i = 0; i < m_edgeGeometry.size(); ++i)
            edgeDists.push_back(m_edgeGeometry[i].signedDistance(p));
        // Create smoothed union geometry around each vertex and then union
        // together
        Real2 dist = 1e5;
        Real2 smoothness = -1.0;
        size_t vtx = 0;
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
#if VERTEX_SMOOTHNESS_MODULATION
                s *= m_vertexSmoothness[u];
#endif
            }
            else {
                jointEdgeDists.clear(), jointEdgeDists.reserve(m_adjEdges[u].size());
                for (size_t ei : m_adjEdges[u])
                    jointEdgeDists.push_back(edgeDists[ei]);
                s = joint->smoothingAmt(p);
#if VERTEX_SMOOTHNESS_MODULATION
                s *= m_vertexSmoothness[u];
#endif
                jdist = SD::exp_smin_reparam_accurate<Real2>(jointEdgeDists, s);
            }
            if (jdist < dist) {
                dist = jdist;
                smoothness = s;
                vtx = u;
            }
        }

        return std::make_pair(smoothness, vtx);
    }

    // Representative cell bounding box (region to be meshed)
    virtual const BBox<Point3D> &boundingBox() const override { return m_bbox; }
    // Choose a different box region to be meshed instead of the default
    // representative cell region (for debugging purposes)
    void setBoundingBox(const BBox<Point3D> &bb) { m_bbox = bb; }

    // Sphere bounding the representative mesh cell, needed for CGAL meshing.
    // Note: CGAL requires the bounding sphere center to lie inside the object.
    virtual void boundingSphere(Point3D &c, double &r) const override {
        auto bbox = boundingBox();
        auto boxCenter = bbox.center();
        // CGAL requires our bounding sphere center to lie within the object, so we
        // find the closest point on the medial axis to the bounding box center.
        c = stripAutoDiff(closestMedialAxisPoint(autodiffCast<Point3<Real>>(boxCenter)));
        // Determine a radius large enough for the bounding box to fit within.
        Point3D boxCorner;
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

    virtual ~PatternSignedDistance() override { }

private:
    ////////////////////////////////////////////////////////////////////////////
    // Signed Distance Evaluation
    ////////////////////////////////////////////////////////////////////////////
    // Additional Real type to support automatic differentiation wrt. p only
    template<typename Real2, bool DebugDerivatives = false>
    Real2 m_signedDistanceImpl(Point3<Real2> p) const {
        p = WMesh::PatternSymmetry::mapToBaseUnit(p);
        std::vector<Real2> edgeDists;
        Real2 closestEdgeDist = 1e5;
        if (m_edgeGeometry.size() == 0) return 1.0;
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
        Real2 dist = combinedJointDistances<Real2, false, DebugDerivatives>(p, edgeDists, closestEdge);

        // TODO: compare distance to each joint against distance to the edge.
        // The joint dist should always be <= edge dist since blending is
        // additive. If both are < edge dist the blending regions overlap and an
        // smin should be used to smooth creases.

        // // Create smoothed union geometry around each vertex and then union
        // // together
        // Real2 dist = 1e5;
        // for (size_t u = 0; u < numVertices(); ++u)
        //     dist = std::min(dist, distToVtxJoint(u));
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

        if (std::isnan(stripAutoDiff(dist)) || std::isinf(stripAutoDiff(dist))) {
            std::cerr << "ERROR, invalid dist: " << dist << "at: " << stripAutoDiff(p).transpose() << std::endl;
            std::cerr << "Edge dists:";
            for (const Real2 &ed : edgeDists) {
                std::cerr << '\t' << ed;
            }
            std::cerr << std::endl;
        }

        assert(!std::isnan(stripAutoDiff(dist)));
        assert(!std::isinf(stripAutoDiff(dist)));

        return dist;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Accelerated Inside/Outside Check
    // For meshers that need only inside-outside queries (e.g. CGAL), we can
    // respond faster if we bypass exact signed distance evaluation.
    ////////////////////////////////////////////////////////////////////////////
    // Running inside-outside tests for the autodiff-typed
    // PatternSignedDistance instantiation is needlessly slow and currently
    // will trigger a static assert in combinedJointDistances--just forbid it.
    template<typename Real2>
    typename std::enable_if<IsAutoDiffType<Real>::value || IsAutoDiffType<Real2>::value, bool>::type
    m_isInsideImpl(Point3<Real2> /* p */) const {
        throw std::runtime_error("Inside/outside check attempted with autodiff class or point.");
    }

    template<typename Real2>
    typename std::enable_if<!(IsAutoDiffType<Real>::value || IsAutoDiffType<Real2>::value), bool>::type
    m_isInsideImpl(Point3<Real2> p) const {
        p = WMesh::PatternSymmetry::mapToBaseUnit(p);
        std::vector<double> edgeDists;
        edgeDists.reserve(m_edgeGeometry.size());
        // Definitely inside if we're inside one of the edges: assumes blending
        // is additive.
        // Note: possibly could be sped up by implementing cheap isInside for
        // edge geometry.
        for (size_t i = 0; i < m_edgeGeometry.size(); ++i) {
            auto dist = stripAutoDiff(m_edgeGeometry[i].signedDistance(autodiffCast<Point3<Real>>(p)));
            if (dist <= 0) return true; // blending is additive
            edgeDists.push_back(dist);
        }

        return combinedJointDistances<Real, true/* accelerated version*/>(p, edgeDists, 0 /* closestEdge */) <= 0;
    }

private:
    const WMesh &m_wireMesh;

    // Bounding box for the meshing cell. Defaults to the representative cell
    // for the symmetry type, but can be changed manually for debugging
    // purposes.
    BBox<Point3D> m_bbox = WMesh::PatternSymmetry::template representativeMeshCell<double>();

    std::vector<typename WMesh::Edge> m_edges; // vertex index pairs for each edge
    std::vector<SD::Primitives::InflatedEdge<Real>> m_edgeGeometry;
    // Quantity between 0 and 1 saying how much to modulate the smoothing of a
    // particular joint (multiplier for m_blendingParams). This is used to
    // prevent bulging of nearly straight joints.
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
