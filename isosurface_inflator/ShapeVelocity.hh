template<typename Real>
class Wireframe {
public:
    // With default parameters
    Wireframe(const std::vector<Point3<Real>> &points,
              const std::vector<std::pair<size_t, size_t>> &edges) {
        m_setGraph(points, edges);
        setParameters(std::vector<Real>(m_points.size(), 0.3),
                      std::vector<Point3<Real>>(m_points.size(), Point3<Real>::Zero()));
    }

    Wireframe(const std::vector<Point3<Real>> &points,
              const std::vector<std::pair<size_t, size_t>> &edges,
              const std::vector<Real> &thicknesses,
              const std::vector<Point3<Real>> &offsets) {
        m_setGraph(points, edges);
        setParameters(thicknesses, offsets);
    }

    void setParameters(const std::vector<Real> &thicknesses,
                       const std::vector<Point3<Real>> &offsets) {
        m_thicknesses = thicknesses;
        m_offsets = offsets;
        m_offsetPoints.clear();
        for (size_t i = 0; i < m_points.size(); ++i)
            m_offsetPoints.emplace_back(m_points[i] + offsets[i]);

        // Vector of connector uses same index as edges
        m_connectorGeometry.clear();
        for (const auto &e : m_edges) {
            m_connectorGeometry.emplace_back(
                    m_offsetPoints[e.first], m_offsetPoints[e.second],
                     m_thicknesses[e.first],  m_thicknesses[e.second]);
        }

        // Compute vertex smoothness:
        // Vertices with intersecting edges are smoothed. This smoothing is
        // ramped up from 0.0 to 1.0 as the minimum incident angle, "theta"
        // shrinks from Pi to Pi/2:
        // smoothness = 0             if theta >= Pi
        //              sin^2(theta)  if Pi/2 < theta < Pi
        //              1.0           if theta < Pi/2
        for (size_t u = 0; u < m_points.size(); ++u) {
            // Min angle over all pairs of connectors
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
                Real angleDeficit1 = uIsP1OfE1 ? m_connectorGeometry.at(e1).angleAtP1()
                                               : m_connectorGeometry.at(e1).angleAtP2();
                // TODO: consider also 2 * angle with symmetry planes...
                // (determine which symmetry planes vertex u lies on)
                for (size_t e2 : m_adjEdges[u]) {
                    if (e2 <= e1) continue;
                    bool uIsP1OfE2 = false;
                    size_t v2 = m_edges[e2].first;
                    if (v2 == u) {
                        uIsP1OfE2 = true;
                        v2 = m_edges[e2].second;
                    }
                    assert(v2 != u);

                    Point3<Real> l1 = m_offsetPoints[v1] - m_offsetPoints[u];
                    Point3<Real> l2 = m_offsetPoints[v2] - m_offsetPoints[u];
                    l1 /= sqrt(l1.squaredNorm()), l2 /= sqrt(l2.squaredNorm());
                    // get angle between edges (in [0, Pi])
                    Real edgeAngle = acos(l1.dot(l2));
                    edgeAngle -= angleDeficit1;
                    edgeAngle -= uIsP1OfE2 ? m_connectorGeometry.at(e2).angleAtP1() : m_connectorGeometry.at(e2).angleAtP2();
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
    }

    void getParameters(std::vector<Real> &thicknesses,
                       std::vector<Point3<Real>> &offsets) const {
        thicknesses = m_thicknesses;
        offsets = m_offsets;
    }

    // Gets the vertex common to edges e1 and e2
    // (or max size_t if none exists)
    size_t commonVertex(size_t e1_i, size_t e2_i) const {
        const auto &e1 = m_edges.at(e1_i);
        const auto &e2 = m_edges.at(e2_i);
        int found = 0;
        size_t v = std::numeric_limits<size_t>::max();
        if ((e1.first == e2.first) || (e1.first == e2.second)) {
            v = e1.first;
            ++found;
        }
        if ((e1.second == e2.first) || (e1.second == e2.second)) {
            v = e1.second;
            ++found;
        }
        if (found == 2) throw std::runtime_error("Invalid edge pair");
        return v;
    }

    size_t numVertices() const { return m_points.size(); }
    size_t valence(size_t v) const { return m_adjEdges.at(v).size(); }

    Real signedDistance(Point3<Real> &p) const {
        std::vector<Real> connectorDists;
        connectorDists.reserve(m_connectorGeometry.size());
        for (const auto &c : m_connectorGeometry)
            connectorDists.push_back(c.signedDistance(p));

        // Create smoothed union geometry around each vertex and then union
        // together
        std::vector<Real> unionedDist;
        std::vector<Real> incidentDists;
        Real dist = 1e5;
        for (size_t u = 0; u < numVertices(); ++u) {
            if (m_adjEdges[u].size() == 0) continue;
            incidentDists.clear();
            for (size_t e : m_adjEdges[u])
                incidentDists.push_back(connectorDists.at(e));
            dist = min(dist, mix(min(incidentDists), exp_smin(incidentDists),
                                 vertexSmoothness(u)));
        }

        return dist;

        // // Optimization: get all connectors whose distances are within the
        // // smoothing radius of each other.
        // // (connector index, distance)
        // std::vector<std::pair<size_t, Real>> dists;
        // for (size_t ci = 0; ci < m_connectorGeometry.size(); ++ci) {
        //     const auto &c = m_connectorGeometry[ci];
        //     Real d = c.signedDistance(p);
        //     if (dists.size()) {
        //         if (d < dists.front()) {
        //             Real newMin = d;
        //             d = dists.front();
        //             dists.front() = newMin;
        //             dists.erase(std::remove_if(dists.begin(), dists.end(),
        //                         [=](Real v) { return (v - newMin) > exp_smin_radius(); }), dists.end());
        //         }
        //         if (d - dists.front() <= exp_smin_radius())
        //             dists.push_back(d);
        //     }
        //     else { dists.push_back(d); }
        //     // dist = std::min(c.signedDistance(p), dist);
        // }
    }

    Real vertexSmoothness(size_t v) const { return m_vertexSmoothness.at(v); }

private:
    void m_setGraph(const std::vector<Point3<Real>> &points,
                    const std::vector<std::pair<size_t, size_t>> &edges) {
        m_points = points;
        // for (auto &p : m_points) p *= 10;

        m_adjEdges.resize(points.size());
        m_adjVertices.resize(points.size());
        for (const auto &e : edges) {
            size_t u = e.first - 1, v = e.second - 1;
            size_t edgeIdx = m_edges.size();
            m_edges.push_back(std::make_pair(u, v));
            m_adjEdges.at(u).push_back(edgeIdx);
            m_adjEdges.at(v).push_back(edgeIdx);
            m_adjVertices.at(u).push_back(v);
            m_adjVertices.at(v).push_back(u);
        }
    }

    void m_determineOrbits() {
    }

    std::vector<std::vector<size_t>>       m_adjEdges;
    std::vector<std::vector<size_t>>       m_adjVertices;
    std::vector<std::pair<size_t, size_t>> m_edges;
    std::vector<Point3<Real>>             m_points;
    std::vector<Point3<Real>>             m_offsetPoints;
    std::vector<Real>                      m_thicknesses;
    std::vector<Point3<Real>>             m_offsets;
    std::vector<Connector<Real>>           m_connectorGeometry;
    // Quantity between 0 and 1 saying how much to smooth intersections between
    // edges incident on a particular vertex
    // Since a given edge pair only has (at most) one vertex in common,
    // specifying this per-vertex makes sense.
    std::vector<Real>                      m_vertexSmoothness;
};

Point3d signedDistanceNormal(const Point3d &p) {
    using adept::adouble;
    typedef Point3<adouble> Pt;
    adept::Stack stack;

    std::vector<Point3<adouble>> pts = { Pt(0.1, 0.1, 0.1), Pt(0.5, 0.5, 0.70), Pt(0.9, 0.9, 0.1) };
    std::vector<adouble> thicknesses = {0.1, 0.25, 0.1};
    std::vector<Point3<adouble>> offsets(pts.size(), Point3<adouble>::Zero());

    std::vector<std::pair<size_t, size_t>> edges = { {1, 2}, {2, 3} }; // 1-indexed
    Wireframe<adouble> pattern(pts, edges);

    pattern.setParameters(thicknesses, offsets);
    Pt x(p.cast<adouble>());
    stack.new_recording();
    for (size_t i = 0; i < 3; ++i) { x[i].set_value(p[i]); }
    adouble dist = pattern.signedDistance(x);
    dist.set_gradient(1);
    stack.reverse();
    Point3<double> distGradX;
    for (size_t i = 0; i < 3; ++i) {
         distGradX[i] =  x[i].get_gradient();
    }
    return distGradX / distGradX.norm();
}

void testShapeVelocity(const vector<Point3d> &points, vector<VectorField<double, 3>> &shapeVelocities) {
    using adept::adouble;
    typedef Point3<adouble> Pt;
    adept::Stack stack;

    std::vector<Point3<adouble>> pts = { Pt(0.1, 0.1, 0.1), Pt(0.5, 0.5, 0.70), Pt(0.9, 0.9, 0.1) };
    std::vector<adouble> thicknesses = {0.1, 0.25, 0.1};
    std::vector<Point3<adouble>> offsets(pts.size(), Point3<adouble>::Zero());

    std::vector<std::pair<size_t, size_t>> edges = { {1, 2}, {2, 3} }; // 1-indexed
    Wireframe<adouble> pattern(pts, edges);

    shapeVelocities.assign(thicknesses.size() + 3 * offsets.size(),
                           VectorField<double, 3>(points.size()));

    size_t pi = 0;
    // NOTE: this could probably be optimized by precomputing the connector
    // objects (not calling pattern.setParameters() every time), but Adept makes
    // this difficult. We'd probably have to manually store the gradients of the
    // computed connector parameters.
    Pt x;
    for (const auto &p : points) {
        stack.new_recording();
        for (size_t i = 0; i < 3; ++i) { x[i].set_value(p[i]); }
        pattern.setParameters(thicknesses, offsets);
        adouble dist = pattern.signedDistance(x);
        dist.set_gradient(1);
        stack.reverse();
        // cout << "Distance: " << dist << endl;
        Point3<double> distGradC1, distGradC2, distGradX;
        
        for (size_t i = 0; i < 3; ++i) {
             distGradX[i] =  x[i].get_gradient();
        }

        std::vector<double> distGradParam;
        distGradParam.reserve(thicknesses.size() + 3 * offsets.size());
        for (const auto &t : thicknesses)
            distGradParam.push_back(t.get_gradient());
        for (const auto &o : offsets) {
            for (size_t c = 0; c < 3; ++c) {
                distGradParam.push_back(o[c].get_gradient());
            }
        }
        for (size_t i = 0; i < distGradParam.size(); ++i)
            shapeVelocities[i](pi) = -distGradParam[i] * distGradX / distGradX.squaredNorm();

        ++pi;
    }
}

