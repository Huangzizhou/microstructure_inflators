////////////////////////////////////////////////////////////////////////////////
// WireMesh.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Pattern graph structure with parameterized embedding.
//      Handles assignment of degrees of freedom to the graph structure based on
//      the requested symmetry and thickness parameter type (both configured by
//      template parameter)
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  06/26/2015 17:54:27
////////////////////////////////////////////////////////////////////////////////
#ifndef WIREMESH_HH
#define WIREMESH_HH

#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>
#include <limits>
#include <set>
#include <queue>

#include <MeshIO.hh>
#include "InflatorTypes.hh"
#include "Symmetry.hh"

// Whether to keep the base cell-adjacent joints in the inflation graph. These
// are needed when their blending regions extend into the base cell.
// Note: slows down meshing/velocity computations.
#define KEEP_BC_ADJ_JOINTS 1

struct TransformedVertex {
    using Point = Point3<double>;
    TransformedVertex(const Point &p, size_t vtx, const Isometry &iso, const Eigen::MatrixXd &pm)
        : pt(p), origVertex(vtx), iso(iso), posMap(pm) { }
    Point pt;
    size_t origVertex;
    Isometry iso;
    Eigen::Matrix3Xd posMap;
    Real sqDistTo(const Point &b) { return (pt - b).squaredNorm(); }
};

struct TransformedEdge {
    using Edge = std::pair<size_t, size_t>;
    TransformedEdge(size_t u, size_t v, size_t oe)
        : e(u, v), origEdge(oe) { }
    UnorderedPair e;
    size_t origEdge;
    // Allow std::map
    bool operator<(const TransformedEdge &b) const { return e < b.e; }
};

enum class ThicknessType { Vertex, Edge };

template<ThicknessType thicknessType_ = ThicknessType::Vertex,
         class Symmetry_ = Symmetry::TriplyPeriodic<>>
class WireMesh {
public:
    typedef Symmetry_ PatternSymmetry;
    typedef Point3<double>            Point;
    typedef std::pair<size_t, size_t> Edge;
    // A vertex adjacent to the base cell: this is represented by the index of
    // its corresponding base cell vertex (first), and the isometry transforming
    // the base cell vertex into it (second).
    typedef std::pair<size_t, Isometry> AdjacentVertex;
    static constexpr ThicknessType thicknessType = thicknessType_;

    WireMesh(const std::string &wirePath) { load(wirePath); }
    WireMesh(const std::vector<MeshIO::IOVertex > &inVertices,
             const std::vector<MeshIO::IOElement> &inElements) {
        set(inVertices, inElements);
    }

    // Embedded graph I/O (OBJ/MSH format)
    void load                  (const std::string &path);
    void save                  (const std::string &path) const;
    void saveBaseUnit          (const std::string &path) const;
    void saveReplicatedBaseUnit(const std::string &path) const;
    void saveInflationGraph    (const std::string &path, std::vector<double> params = std::vector<double>()) const;

    void set(const std::vector<MeshIO::IOVertex > &inVertices,
             const std::vector<MeshIO::IOElement> &inElements);

    size_t numVertices    () const { return m_fullVertices.size(); }
    size_t numEdges       () const { return m_fullEdges   .size(); }
    size_t numBaseVertices() const { return m_baseVertices.size(); }
    size_t numBaseEdges   () const { return m_baseEdges   .size(); }

    // There is a thickness parameter for each base vertex or base edge
    // depending on thicknessType.
    size_t numThicknessParams() const {
        switch (thicknessType) {
            case ThicknessType::Edge:   return m_baseEdges.size();
            case ThicknessType::Vertex: return m_baseVertices.size();
            default: assert(false);
        }
    }
    // Positional parameters always live on the base vertices, and they are
    // determined by the base vertex positioner.
    size_t numPositionParams() const {
        size_t posParams = 0;
        for (const auto &bvp : m_baseVertexPositioners)
            posParams += bvp.numDoFs();
        return posParams;
    }

    // There is a single blending parameter per base vertex.
    size_t numBlendingParameters() const { return numBaseVertices(); }

    size_t numParams() const { return numThicknessParams() + numPositionParams() + numBlendingParameters(); }

    // Determine the position parameters from the original embedded graph
    // positions.
    std::vector<double> defaultPositionParams() const {
        size_t np = numPositionParams();
        std::vector<double> positionParams(np);
        size_t paramOffset = 0;
        for (size_t i = 0; i < m_baseVertices.size(); ++i) {
            m_baseVertexPositioners[i].getDoFsForPoint(m_baseVertices[i],
                    &positionParams[paramOffset]);
            paramOffset += m_baseVertexPositioners[i].numDoFs();
        }
        assert(paramOffset == np);
        return positionParams;
    }

    // TODO: MAKE CONFIGURABLE.
    std::vector<double> defaultThicknessParams() const {
        return std::vector<double>(numThicknessParams(), 0.07);
    }

    std::vector<double> defaultBlendingParams() const {
        return std::vector<double>(numBlendingParameters(), 0.01);
    }

    // Position parameters come first, followed by thickness and blending
    bool isPositionParam(size_t p) const {
        if (p >= numParams()) throw std::runtime_error("Invalid parameter index.");
        return p < numPositionParams();
    }
    bool isThicknessParam(size_t p) const {
        if (p >= numParams()) throw std::runtime_error("Invalid parameter index.");
        return (p >= numPositionParams()) && (p < numPositionParams() + numThicknessParams());
    }
    bool isBlendingParam(size_t p) const {
        if (p >= numParams()) throw std::runtime_error("Invalid parameter index.");
        return p >= numPositionParams() + numThicknessParams();
    }

    // Position parameters come first, followed by thickness and blending
    std::vector<double> defaultParameters() const {
        std::vector<double> params;
        params.reserve(numParams());
        params = defaultPositionParams();
        auto dtp = defaultThicknessParams(); params.insert(params.end(), dtp.begin(), dtp.end());
        auto dbp =  defaultBlendingParams(); params.insert(params.end(), dbp.begin(), dbp.end());
        return params;
    }

    // The inflation graph includes all vertices and edges in the base symmetry
    // unit plus adjacent edges and vertices. The adjacent subgraph is needed to
    // generate the correct inflated joint geometry at the interface nodes.
    // The "points" returned are after the positional parameters have been
    // applied, and the thickness parameters are decoded from "params" into the
    // "thicknesses" vector for convenience: these will correspond to either
    // entries in "points" or "edges" depending on thicknessType.
    // The blending parameters are also decoded into the per-vertex
    // "blendingParams" vector.
    template<typename Real>
    void inflationGraph(const std::vector<Real> &params,
                        std::vector<Point3<Real>> &points,
                        std::vector<Edge> &edges,
                        std::vector<Real> &thicknesses,
                        std::vector<Real> &blendingParams) const {
        if (params.size() != numParams())
            throw std::runtime_error("Invalid number of params.");

        edges = m_inflEdge;

        // Set the base graph positions using params
        std::vector<Point3<Real>> baseGraphPos;
        baseGraphPos.reserve(m_baseVertices.size());
        for (size_t i = 0, pOffset = 0; i < m_baseVertices.size(); ++i) {
            const auto &pos = m_baseVertexPositioners[i];
            baseGraphPos.push_back(pos.template getPosition<Real>(&params[pOffset]));
            pOffset += pos.numDoFs();
        }

        points.clear(), points.reserve(m_inflVtx.size());
        // Determine the inflation graph positions from the base positions.
        for (const TransformedVertex &iv : m_inflVtx)
            points.push_back(iv.iso.apply(baseGraphPos.at(iv.origVertex)));

        // Determine thicknesses and blending parameters
        if (thicknessType == ThicknessType::Edge)
            throw std::runtime_error("Edge thickness currently unimplemented"); // will be easy to implement with m_inflEdgeOrigin

        thicknesses   .clear(), thicknesses   .reserve(m_inflVtx.size());
        blendingParams.clear(), blendingParams.reserve(m_inflVtx.size());

        const size_t thicknessVarOffsets = numPositionParams();
        const size_t blendingVarOffsets = thicknessVarOffsets + numThicknessParams();
        for (const TransformedVertex &iv : m_inflVtx) {
            thicknesses   .push_back(params.at(thicknessVarOffsets + iv.origVertex));
            blendingParams.push_back(params.at( blendingVarOffsets + iv.origVertex));
        }
    }

    // Construct (stitched) replicated graph along with the maps from parameters
    // to vertex positions/thicknesses/blending parameters
    // Note: these maps operate on "homogenous parameters" (i.e. the vector of
    // parameters with 1 appended).
    // TODO: make sure we have all parameters set up before this is called from
    // the constructor (to create inflation graph)
    void replicatedGraph(const std::vector<Isometry> &isometries,
                         std::vector<TransformedVertex> &outVertices,
                         std::vector<TransformedEdge  > &outEdges) {
        std::vector<TransformedVertex> rVertices;
        std::set<TransformedEdge>      rEdges;

        // Determine the index of the first position parameter for each vertex
        std::vector<size_t> vtxPosParamOffset;
        size_t posParams = 0;
        const size_t nparams = numParams();
        for (const auto &bvp : m_baseVertexPositioners) {
            vtxPosParamOffset.push_back(posParams);
            posParams += bvp.numDoFs();
        }

        for (size_t ei = 0; ei < m_baseEdges.size(); ++ei) {
            const auto &e = m_baseEdges[ei];
            const size_t u = e.first, v = e.second;

            auto pu = m_baseVertices.at(u),
                 pv = m_baseVertices.at(v);

            Eigen::Matrix3Xd uPosMap = m_baseVertexPositioners.at(u).getPositionMap(nparams, vtxPosParamOffset.at(u));
            Eigen::Matrix3Xd vPosMap = m_baseVertexPositioners.at(v).getPositionMap(nparams, vtxPosParamOffset.at(v));

            for (const auto &isometry : isometries) {
                auto mappedPu = isometry.apply(pu);
                auto mappedPv = isometry.apply(pv);

                const bool hasTranslation = isometry.hasTranslation();

                size_t mappedUIdx = std::numeric_limits<size_t>::max(),
                       mappedVIdx = std::numeric_limits<size_t>::max();
                // Search for repeated vertices.
                for (size_t i = 0; i < rVertices.size(); ++i) {
                    if (rVertices[i].sqDistTo(mappedPu) < PatternSymmetry::tolerance) {
                        if (!hasTranslation) // Translations will give different position maps, but this should be fine for orthotropic patterns
                            assert((rVertices[i].posMap - isometry.xformMap(uPosMap)).squaredNorm() < 1e-8);
                        assert(rVertices[i].origVertex == u); // Note: will fail on triply periodic; need constraints
                        mappedUIdx = i;
                    }
                    if (rVertices[i].sqDistTo(mappedPv) < PatternSymmetry::tolerance) {
                        if (!hasTranslation)
                            assert((rVertices[i].posMap - isometry.xformMap(vPosMap)).squaredNorm() < 1e-8);
                        assert(rVertices[i].origVertex == v); // Note: will fail on triply periodic; need constraints
                        mappedVIdx = i;
                    }
                }
                // Add the mapped vertices if they are distinct from the
                // existing ones.
                if (mappedUIdx > rVertices.size()) {
                    mappedUIdx = rVertices.size();
                    rVertices.emplace_back(mappedPu, u, isometry, isometry.xformMap(uPosMap));
                }
                if (mappedVIdx > rVertices.size()) {
                    mappedVIdx = rVertices.size();
                    rVertices.emplace_back(mappedPv, v, isometry, isometry.xformMap(vPosMap));
                }

                TransformedEdge xfEdge(mappedUIdx, mappedVIdx, ei);
                auto ret = rEdges.insert(xfEdge);
                if (ret.second == false) // reflected already exists; verify it's from the same base edge
                    assert(ret.first->origEdge == ei); // Note: will fail on triply periodic
            }
        }

        outVertices = rVertices;
        outEdges.clear(), outEdges.reserve(rEdges.size());
        for (const auto &re : rEdges)
            outEdges.push_back(re);
    }

    // Extract the replicated graph needed to validate printability (i.e. self
    // supporting check) and to determine the self supporting constraints for
    // the optimizer.
    // For general periodic patterns, this is the full cell.
    // For cubic/orthotropic patterns, this is one full "column" of the period
    // cell (x, y in positive quadrant, z in [-1, 1]).
    //
    // Assumes ThicknessType::Vertex for now.
    // The thickness and position of each vertex is decoded from the "params"
    // vector into "points" and "thickness vars" for the printability check.
    //
    // For determining the self-supporting constraints, the linear map from
    // pattern parameters to vertex positions is encoded in the "positionMap" and
    // "thicknessMap" matrices.
    // These matrices will actually be sparse, but using a dense representation
    // shouldn't be a bottleneck.
    // The last column of these maps is a constant offset (like in homogeneous
    // coordinates)
    void printabilityGraph(const std::vector<Real> &params,
                           std::vector<Point3<Real>> &points,
                           std::vector<Edge> &edges,
                           std::vector<size_t> &thicknessVars,
                           std::vector<Eigen::Matrix3Xd> &positionMaps)
    {
        if (params.size() != numParams())
            throw std::runtime_error("Invalid number of params.");

        // Reflect base graph into the "printability column." This is done by
        // applying all permutation isometries, but no translation and only the
        // z-axis reflection
        std::vector<Isometry> printabiltyIsometries;
        for (const auto &iso : PatternSymmetry::symmetryGroup()) {
            if (iso.hasTranslation() ||
                iso.hasReflection(Symmetry::Axis::X) ||
                iso.hasReflection(Symmetry::Axis::Y)) {
                continue;
            }
            printabiltyIsometries.push_back(iso);
        }

        std::vector<TransformedVertex> rVertices;
        std::vector<TransformedEdge>   rEdges;
        replicatedGraph(printabiltyIsometries,
                        rVertices, rEdges);

        points       .clear(); points       .reserve(rVertices.size());
        edges        .clear(); edges        .reserve(   rEdges.size());
        thicknessVars.clear(); thicknessVars.reserve(rVertices.size());
        positionMaps .clear(); positionMaps .reserve(rVertices.size());

        // Use position map to place points; we need to construct homogenous
        // param vector.
        Eigen::VectorXd paramVec(params.size() + 1);
        for (size_t i = 0; i < params.size(); ++i) paramVec[i] = params[i];
        paramVec[params.size()] = 1.0;

        static_assert(thicknessType_ == ThicknessType::Vertex, "Only Per-Vertex Thickness Supported");
        const size_t tvOffset = numPositionParams();
        for (const auto &rv : rVertices) {
            points.emplace_back(rv.posMap * paramVec);
            thicknessVars.emplace_back(tvOffset + rv.origVertex);
            positionMaps.emplace_back(rv.posMap);
        }

        for (const auto &re : rEdges)
            edges.push_back({re.e[0], re.e[1]});
    }

    bool isPrintable(const std::vector<Real> &params) const {
        if (thicknessType_ != ThicknessType::Vertex)
            throw std::runtime_error("Only per-vertex thickness printability implemented currently.");

        std::vector<Point3<Real>> points;
        std::vector<Edge> edges;
        std::vector<size_t> thicknessVars;
        std::vector<Eigen::Matrix3Xd> positionMaps;
        printabilityGraph(params, points, edges, thicknessVars, positionMaps);

        assert(thicknessVars.size() == points.size());

        std::vector<Real> thicknesses;
        for (size_t tv : thicknessVars)
            thicknesses.push_back(params.at(tv));

        std::vector<Real> zCoords;
        for (const auto &pt : points)
            zCoords.push_back(pt[2]);

        // Build adjacency list
        std::vector<std::vector<size_t>> adj(points.size());
        for (const auto &e : edges) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }

        std::vector<bool> supported(points.size(), false);
        double tol = PatternSymmetry::tolerance;

        // Should actually be the orthotropic base cell bottom
        Real minZ = *std::min_element(zCoords.begin(), zCoords.end());

        // Subtract radius from each vertex to get height of vertex sphere's
        // bottom
        for (size_t i = 0; i < points.size(); ++i)
            zCoords[i] -= thicknesses[i];

        std::queue<size_t> bfsQueue;
        for (size_t i = 0; i < zCoords.size(); ++i) {
            if (zCoords[i] - minZ < tol) {
                bfsQueue.push(i);
                supported[i] = true;
            }
        }

        while (!bfsQueue.empty()) {
            size_t u = bfsQueue.front();
            bfsQueue.pop();
            for (size_t v : adj[u]) {
                if (!supported[v] && (zCoords[v] - zCoords[u] >= -tol)) {
                    supported[v] = true;
                    bfsQueue.push(v);
                }
            }
        }

        bool result = true;
        for (bool s : supported) result &= s;
        return result;
    }

    // The self-supporting printability constraints for per-vertex thickness
    // patterns take the form of inequality constraints on the position and
    // thickness parameters. These inequalities are linear apart from a min()
    // operation on the supporting candidates' z coordinates (all vertices are
    // supported if each is above the minimum of its neighbors.)
    // This min() is applied based on "params" so that a set of linear
    // inequality constraints is returned.
    // TODO: look up if these are actually "affine" constraints; is that a
    // thing?
    // The constraints are in the form:
    //      C [p] >= 0
    //        [1]
    // where p is the parameter vector and C is the constraint matrix returned.
    //
    // Perfectly horizontal features are forbidden because they make
    // printability a global condition. Thus each vertex is required to be
    // at least "epsilon" above its neighbors.
    // Note: this should always be possible in the per-vertex-thickness model,
    // but the per-edge-thickness model does require truly horizontal edges in
    // the orthotropic symmetry.
    Eigen::MatrixXd selfSupportingConstraints(
            const std::vector<Real> &params,
            Real epsilon) const
    {
        if (thicknessType_ != ThicknessType::Vertex)
            throw std::runtime_error("Only per-vertex thickness printability implemented currently.");

        std::vector<Point3<Real>> points;
        std::vector<Edge> edges;
        std::vector<size_t> thicknessVars;
        std::vector<Eigen::Matrix3Xd> positionMaps;
        printabilityGraph(params, points, edges, thicknessVars, positionMaps);

        assert(positionMaps.size() == thicknessVars.size());

        std::vector<Real> zCoords;
        for (const auto &pt : points)
            zCoords.push_back(pt[2]);

        Real minZ = *std::min_element(zCoords.begin(), zCoords.end());

        // Subtract the radius from each vertex position to get the sphere
        // bottom point; now positionMaps map pattern parameters to the sphere's
        // lowest point.
        for (size_t i = 0; i < positionMaps.size(); ++i) {
            positionMaps[i](2, thicknessVars[i]) -= 1.0;
            zCoords[i] -= params.at(thicknessVars[i]);
        }

        // Build adjacency list
        std::vector<std::vector<size_t>> adj(points.size());
        for (const auto &e : edges) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }

        double tol = PatternSymmetry::tolerance;
        // Create a constraint for every vertex
        std::vector<Eigen::VectorXd> constraints;
        for (size_t u = 0; u < points.size(); ++u) {
            // ... except those on the build platform
            if (zCoords[u] - minZ < tol)
                continue;

            // Find lowest neighbor
            Real minHeight = safe_numeric_limits<  Real>::max();
            size_t lowestNeighbor  = safe_numeric_limits<size_t>::max();
            for (size_t v : adj[u]) {
                if (zCoords[v] < minHeight) {
                    minHeight = zCoords[v];
                    lowestNeighbor = v;
                }
            }

            // Don't allow perfectly horizontal features since these require
            // global printability conditions.
            // u - lowestNeighbor           >= epsilon  <==>
            // u - lowestNeighbor - epsilon >= 0
            constraints.push_back(positionMaps[u].row(2) -
                                  positionMaps[lowestNeighbor].row(2));
            assert(constraints.back().cols() == params.size() + 1);
            constraints.back()[params.size()] -= epsilon;
        }

        // Detect and remove "constant" constraints not acting on variables;
        // these must be satisfied
        std::vector<Eigen::VectorXd> prunedConstraints;
        for (const auto &c : constraints) {
            bool hasVar = false;
            for (size_t p = 0; p < params.size(); ++p) {
                if (std::abs(c[p]) > tol) {
                    hasVar = true;
                    break;
                }
            }
            if (hasVar) {
                prunedConstraints.emplace_back(c);
            }
            else {
                // Constraints c [p] >= 0 not acting on vars must be satisfied.
                //               [1]
                if (c[params.size()] < -tol)
                    throw std::runtime_error("Infeasible: constant constraint unsatisfied");
            }
        }

        // TODO: think about reducing the constraint system

        Eigen::MatrixXd C(prunedConstraints.size(), params.size() + 1);
        for (size_t i = 0; i < prunedConstraints.size(); ++i)
            C.row(i) = prunedConstraints[i];
        return C;
    }

private:

    // All vertex/edges of the pattern
    std::vector<Point>  m_fullVertices;
    std::vector<Edge>   m_fullEdges;
    // The distinct vertices/edges of the pattern modulo symmetry
    std::vector<Point>  m_baseVertices;
    std::vector<Edge>   m_baseEdges;

    std::vector<decltype(PatternSymmetry::nodePositioner(Point()))> m_baseVertexPositioners;

    // The inflation graph, with each vertex/edge linked back to the
    // corresponding originating vertex/edge in the base cell.
    // This graph consists of only the edges needed to properly define the
    // geometry inside the base cell (i.e. all base cell edges, edges incident
    // the base cell. For cases where blending regions for joints outside the
    // base cell can extend inside, we also include edges adjacent to the
    // base-cell-incident edges).
    std::vector<Edge>              m_inflEdge;
    std::vector<size_t>            m_inflEdgeOrigin;
    std::vector<TransformedVertex> m_inflVtx;

    // Find the base vertex within symmetry tolerance of p
    // (or throw an exception if none exists).
    size_t m_findBaseVertex(const Point &p) const;
};

#include "WireMesh.inl"

#endif /* end of include guard: WIREMESH_HH */
