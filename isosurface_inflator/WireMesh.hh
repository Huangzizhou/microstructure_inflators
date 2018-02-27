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
#include <array>
#include <algorithm>
#include <numeric>

#include <MeshIO.hh>
#include <Utilities/apply.hh>
#include "InflatorTypes.hh"
#include "Symmetry.hh"
#include "AutomaticDifferentiation.hh"

struct TransformedVertex {
    using Point = Point3<double>;
    TransformedVertex(const Point &p, size_t vtx, const Isometry &iso, const Eigen::MatrixXd &pm)
            : pt(p), origVertex(vtx), iso(iso), posMap(pm) { }
    Point pt;
    size_t origVertex;
    Isometry iso; // map from origin vertex to this reflected copy
    Eigen::Matrix3Xd posMap; // matrix rep. of map from params to vtx position
                             // (map is affine; last col is a const translation)
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

// Stores the positions of all of a base vertex's variables in the
// full parameter vector.
// Variables can be vector-valued (e.g., position vars), in which case this
// struct stores only the position of the first component--the remaining
// components follow contiguously.
// Note that base vertices can share variables (e.g., for patterns with only
// triply periodic symmetry).
struct BaseVtxVarOffsets { size_t position, thickness, blending; bool ignored; };

// Parts of WireMesh's interface that can be implemented without the
// PatternSymmetry template parameter (or nan be made virtual)
class WireMeshBase {
public:
    using Point = Point3<double>;
    using Edge = std::pair<size_t, size_t>;

    virtual void set(const std::vector<MeshIO::IOVertex > &inVertices,
                     const std::vector<MeshIO::IOElement> &inElements) = 0;

    // Embedded graph I/O (OBJ/MSH format)
    void load                  (const std::string &path);
    void save                  (const std::string &path) const;
    void saveBaseUnit          (const std::string &path) const;
    void saveReplicatedBaseUnit(const std::string &path) const;
    void savePeriodCellGraph   (const std::string &path) const;
    void saveInflationGraph    (const std::string &path, std::vector<double> params = std::vector<double>()) const;

    size_t numVertices    () const { return m_fullVertices.size(); }
    size_t numEdges       () const { return m_fullEdges   .size(); }
    size_t numBaseVertices() const { return m_baseVertices.size(); }
    size_t numBaseEdges   () const { return m_baseEdges   .size(); }

    // Positional parameters always live on the base vertices, and they are
    // determined by the base vertex positioner.
    size_t numPositionParams() const { return m_numPositionParams; }

    ThicknessType thicknessType() const { return m_thicknessType; }

    // There is a thickness parameter for each base vertex or base edge
    // depending on thicknessType.
    size_t numThicknessParams() const {
        switch (thicknessType()) {
            case ThicknessType::Edge:   return numBaseEdges();
            case ThicknessType::Vertex: return m_numIndepBaseVertices;
            default: assert(false);
        }
    }

    // There is a single blending parameter per base vertex.
    size_t numBlendingParameters() const { return m_numIndepBaseVertices; }

    size_t numParams() const { return numPositionParams() + numThicknessParams() + numBlendingParameters(); }
    size_t numFilteredInParams() const {
        if (!m_useFilterMask)
            return numParams();

        return std::count(m_parametersMask.begin(), m_parametersMask.begin()+numParams(), false);
    }

    virtual std::vector<double> defaultPositionParams() const = 0;

    // Position parameters come first, followed by thickness and blending
    std::vector<double> defaultParameters(double thickness = 0.07) const {
        std::vector<double> params;
        params.reserve(numParams());
        params = defaultPositionParams();
        auto dtp = defaultThicknessParams(thickness); params.insert(params.end(), dtp.begin(), dtp.end());
        auto dbp =  defaultBlendingParams(); params.insert(params.end(), dbp.begin(), dbp.end());
        return params;
    }

    // TODO: test that parametersForPeriodCellGraph and periodCellGraph are inverses
    virtual void parametersForPeriodCellGraph(
        const std::vector<Point3<double>> &points,
        const std::vector<Edge> &edges,
        const std::vector<double> &thicknesses,
        const std::vector<double> &blendingParams,
        std::vector<double> &params) const = 0;

    // TODO: MAKE CONFIGURABLE.
    std::vector<double> defaultThicknessParams(double thickness = 0.07) const {
        return std::vector<double>(numThicknessParams(), thickness);
    }

    std::vector<double> defaultBlendingParams() const {
        return std::vector<double>(numBlendingParameters(), 0.01);
    }

    // Position parameters come first, followed by thickness and blending
    void validateParamIdx(size_t p) const {
        if (p >= numParams())
            throw std::runtime_error("Invalid parameter index"); }

    bool  isPositionParam(size_t p) const {
        if (m_useFilterMask) {
            p = m_filteredParamsToOriginal[p];
        }

        validateParamIdx(p);

        return p < numPositionParams();
    }
    bool isThicknessParam(size_t p) const {
        if (m_useFilterMask) {
            p = m_filteredParamsToOriginal[p];
        }

        validateParamIdx(p);
        return (p >= numPositionParams()) && (p < numPositionParams() + numThicknessParams());
    }
    bool  isBlendingParam(size_t p) const {
        if (m_useFilterMask) {
            p = m_filteredParamsToOriginal[p];
        }

        validateParamIdx(p);
        return p >= numPositionParams() + numThicknessParams();
    };

    virtual bool isPrintable(const std::vector<Real> &params, bool verticalInterfacesSupported = false) const = 0;

    virtual void inflationGraph(const std::vector<double> &params,
                                std::vector<Point3<double>> &points,
                                std::vector<Edge> &edges,
                                std::vector<double> &thicknesses,
                                std::vector<double> &blendingParams) const = 0;

    // Full period cell graph in its default configuration.
    void periodCellGraph(std::vector<Point> &points,
                         std::vector<Edge>  &edges) const {
        edges  = m_periodCellEdge;
        points = apply(m_periodCellVtx, [](const TransformedVertex &tv) { return tv.pt; });
    }

    // Full period cell graph and associated thickness/blending params under
    // parameter values "params".
    virtual void periodCellGraph(const std::vector<double> &params,
                         std::vector<Point>  &points,
                         std::vector<Edge>   &edges,
                         std::vector<double> &thicknesses,
                         std::vector<double> &blendingParams) const = 0;

    // Construct (stitched) replicated graph along with the maps from parameters
    // to vertex positions/thicknesses/blending parameters
    // Note: these maps operate on "homogeneous parameters" (i.e. the vector of
    // parameters with 1 appended).
    virtual void replicatedGraph(const std::vector<Isometry> &isometries,
                         std::vector<TransformedVertex> &outVertices,
                         std::vector<TransformedEdge  > &outEdges) const = 0;

    virtual std::vector<Isometry> symmetryGroup() const = 0;

protected:
    // All vertex/edges of the pattern
    std::vector<Point>  m_fullVertices;
    std::vector<Edge>   m_fullEdges;
    // The distinct vertices/edges of the pattern modulo symmetry
    std::vector<Point>  m_baseVertices;
    std::vector<Edge>   m_baseEdges;

    // For certain pattern symmetries (e.g. TriplyPeriodic), some vertices in
    // the base graph must be constrained to have the same
    // position/thickness/radius variables as others. We do this by declaring
    // certain vertices to be "independent" and others to be "dependent."
    // Dependent vertices share variables with their single paired independent
    // vertex (m_indepVtxForBaseVtx[*]).
    size_t m_numIndepBaseVertices,
           m_numDepBaseVertices;
    std::vector<size_t> m_indepVtxForBaseVtx;

    size_t m_numPositionParams;
    // Index into the parameter vector of each base vertex's parameters.
    std::vector<BaseVtxVarOffsets> m_baseVertexVarOffsets;

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

    // The full period cell graph
    std::vector<Edge>              m_periodCellEdge;
    std::vector<TransformedVertex> m_periodCellVtx;

    // The printability graph
    std::vector<Edge>              m_printGraphEdge;
    std::vector<TransformedVertex> m_printGraphVtx;

    ThicknessType m_thicknessType = ThicknessType::Vertex;

    virtual double m_tolerance() const = 0;

    std::vector<bool> m_parametersMask;
    std::vector<Real> m_originalParameters;
    std::vector<int> m_filteredParamsToOriginal;
    bool m_useFilterMask = false;
};

template<class Symmetry_ = Symmetry::TriplyPeriodic<>>
class WireMesh : public WireMeshBase {
public:
    using PatternSymmetry = Symmetry_;

    std::vector<int> generateFilterMap(const std::vector<bool> &parametersMask) {
        int filteredIdx = 0;
        int originalIdx = 0;
        std::vector<int> result;

        for (; originalIdx < parametersMask.size(); originalIdx++) {
            if (!parametersMask[originalIdx]) {
                result.push_back(originalIdx);
            }
        }

        return result;
    }

    WireMesh(const std::string &wirePath, size_t inflationNeighborhoodEdgeDist = 2)
        : m_inflationNeighborhoodEdgeDist(inflationNeighborhoodEdgeDist) { load(wirePath); }
    WireMesh(const std::string &wirePath, const std::vector<bool> &parametersMask, const std::vector<double> originalParameters, size_t inflationNeighborhoodEdgeDist = 2)
            : m_inflationNeighborhoodEdgeDist(inflationNeighborhoodEdgeDist) {
        m_parametersMask = parametersMask;
        m_originalParameters = originalParameters;
        m_filteredParamsToOriginal = generateFilterMap(parametersMask);
        m_useFilterMask = true; load(wirePath);
    }
    WireMesh(const std::vector<MeshIO::IOVertex > &inVertices,
             const std::vector<MeshIO::IOElement> &inElements,
             size_t inflationNeighborhoodEdgeDist = 2)
        : m_inflationNeighborhoodEdgeDist(inflationNeighborhoodEdgeDist) { set(inVertices, inElements); }

    virtual std::vector<Isometry> symmetryGroup() const override {
        return PatternSymmetry::symmetryGroup();
    }

    virtual void set(const std::vector<MeshIO::IOVertex > &inVertices,
                     const std::vector<MeshIO::IOElement> &inElements) override;

    // Determine the position parameters from the original embedded graph
    // positions.
    virtual std::vector<double> defaultPositionParams() const override {
        std::vector<double> positionParams(numPositionParams());
        for (size_t i = 0; i < m_baseVertices.size(); ++i) {
            if (m_baseVertexPositioners[i].numDoFs() == 0) continue;
            size_t offset = m_baseVertexVarOffsets[i].position;
            assert(offset < positionParams.size());
            m_baseVertexPositioners[i].getDoFsForPoint(m_baseVertices[i],
                                                       &positionParams[offset]);
        }

        return positionParams;
    }

    virtual void parametersForPeriodCellGraph(
            const std::vector<Point3<double>> &points,
            const std::vector<Edge> &/* edges */,
            const std::vector<double> &thicknesses,
            const std::vector<double> &blendingParams,
            std::vector<double> &params) const override {
        params.assign(numParams(), 0);
        assert(points.size() == m_periodCellVtx.size());
        for (size_t i = 0; i < m_periodCellVtx.size(); ++i) {
            const auto &pcv = m_periodCellVtx[i];
            if (!pcv.iso.isIdentity()) continue;
            size_t bi = pcv.origVertex;
            size_t offset = m_baseVertexVarOffsets[bi].position;
            assert(offset < params.size());
            m_baseVertexPositioners[bi].getDoFsForPoint(points.at(i),
                                                        &params[offset]);
            params.at(m_baseVertexVarOffsets[bi].thickness) = thicknesses.at(i);
            params.at(m_baseVertexVarOffsets[bi].blending) = blendingParams.at(i);
        }
    }

    virtual void inflationGraph(const std::vector<double> &params,
                                std::vector<Point3<double>> &points,
                                std::vector<Edge> &edges,
                                std::vector<double> &thicknesses,
                                std::vector<double> &blendingParams) const override {
        inflationGraph<double>(params, points, edges, thicknesses, blendingParams);
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
        m_getGraphForParameters(m_inflVtx, m_inflEdge, params, points, edges,
                                thicknesses, blendingParams);
    }

    // Full period cell graph and associated thickness/blending params under
    // parameter values "params".
    virtual void periodCellGraph(const std::vector<double> &params,
                         std::vector<Point3<double>> &points,
                         std::vector<Edge>   &edges,
                         std::vector<double> &thicknesses,
                         std::vector<double> &blendingParams) const override {
        m_getGraphForParameters(m_periodCellVtx, m_periodCellEdge, params, points, edges,
                                thicknesses, blendingParams);
    }

    // Construct (stitched) replicated graph along with the maps from parameters
    // to vertex positions/thicknesses/blending parameters
    // Note: these maps operate on "homogeneous parameters" (i.e. the vector of
    // parameters with 1 appended).
    virtual void replicatedGraph(const std::vector<Isometry> &isometries,
                         std::vector<TransformedVertex> &outVertices,
                         std::vector<TransformedEdge  > &outEdges) const override;

    // Extract the replicated graph needed to validate printability (i.e. self
    // supporting check) and to determine the self supporting constraints for
    // the optimizer.
    // For general periodic patterns, this is the full cell.
    // For cubic/orthotropic patterns, this is one full "column" of the period
    // cell (x, y in positive quadrant, z in [-1, 1]).
    //
    // Assumes ThicknessType::Vertex for now.
    //
    // For determining the self-supporting constraints, the linear map from
    // pattern parameters to vertex positions is encoded in the "positionMap" and
    // "thicknessMap" matrices.
    // These matrices will be sparse, but using a dense representation
    // shouldn't be a bottleneck.
    // The last column of these maps is a constant offset (as in homogeneous
    // coordinates)
    void printabilityGraph(std::vector<Edge> &edges,
                           std::vector<size_t> &thicknessVars,
                           std::vector<Eigen::Matrix3Xd> &positionMaps) const
    {
        if (thicknessType() != ThicknessType::Vertex)
            throw std::runtime_error("Only Per-Vertex Thickness Supported");

        // Decode from cached printability graph.
        thicknessVars.clear(); thicknessVars.reserve(m_printGraphVtx.size());
        positionMaps .clear(); positionMaps .reserve(m_printGraphVtx.size());
        for (const auto &pv : m_printGraphVtx) {
            thicknessVars.emplace_back(m_baseVertexVarOffsets.at(pv.origVertex).thickness);
            positionMaps.emplace_back(pv.posMap);
        }

        edges = m_printGraphEdge;
    }

    // Same as above, but also get the position of each vertex according to
    // "params"
    void printabilityGraph(const std::vector<Real> &params,
                           std::vector<Point3<Real>> &points,
                           std::vector<Edge> &edges,
                           std::vector<size_t> &thicknessVars,
                           std::vector<Eigen::Matrix3Xd> &positionMaps) const
    {
        if (params.size() != numParams())
            throw std::runtime_error("Invalid number of params.");

        printabilityGraph(edges, thicknessVars, positionMaps);

        // Use position map to place points; we need to construct homogeneous
        // param vector.
        Eigen::VectorXd paramVec(params.size() + 1);
        for (size_t i = 0; i < params.size(); ++i) paramVec[i] = params[i];
        paramVec[params.size()] = 1.0;

        points.clear(); points.reserve(positionMaps.size());
        for (const auto &pm : positionMaps)
            points.emplace_back(pm * paramVec);
    }

    // verticalInterfacesSupported: whether to assume the vertical interface nodes are supported
    // When working with tilings (StitchedWireMesh), the averaging stitching operation
    // break printability. However, this operation never makes the interface nodes
    // unprintable (since the lowest interface node was supported before
    // averaging and only moves up during averaging). To correctly analyze printability
    // in this case (while looking only at a single period cell at a time), we
    // must manually specify that the vertical interfaces are supported.
    bool isPrintable(const std::vector<Real> &params, bool verticalInterfacesSupported = false) const override {
        if (thicknessType() != ThicknessType::Vertex)
            throw std::runtime_error("Only Per-Vertex Thickness Supported");

        std::vector<Edge> edges;
        std::vector<Point> points;
        std::vector<size_t> thicknessVars;
        std::vector<Eigen::Matrix3Xd> positionMaps;
        printabilityGraph(params, points, edges, thicknessVars, positionMaps);
        const size_t numPoints = points.size();

        // Get the heights of the vertex spheres' bottoms
        std::vector<Real> zCoords; zCoords.reserve(numPoints);
        for (size_t i = 0; i < numPoints; ++i)
            zCoords.push_back(points[i][2] - params.at(thicknessVars[i]));

        // Build adjacency list
        std::vector<std::vector<size_t>> adj(numPoints);
        for (const auto &e : edges) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }

        std::vector<bool> supported(numPoints, false);
        const double tol = m_tolerance();

        // Should actually be the base cell bottom
        Real minZ = *std::min_element(zCoords.begin(), zCoords.end());

        std::queue<size_t> bfsQueue;
        for (size_t i = 0; i < numPoints; ++i) {
            bool s = zCoords[i] < minZ + tol;
            s |= verticalInterfacesSupported &&
                 ((std::abs(points[i][0] - 1) < PatternSymmetry::tolerance) ||
                  (std::abs(points[i][0] + 1) < PatternSymmetry::tolerance) ||
                  (std::abs(points[i][1] - 1) < PatternSymmetry::tolerance) ||
                  (std::abs(points[i][1] + 1) < PatternSymmetry::tolerance));
            if (s) {
                bfsQueue.push(i);
                supported[i] = true;
            }
        }

        while (!bfsQueue.empty()) {
            size_t u = bfsQueue.front();
            bfsQueue.pop();
            for (size_t v : adj[u]) {
                if (!supported[v] && (zCoords[v] >= zCoords[u] - tol)) {
                    supported[v] = true;
                    bfsQueue.push(v);
                }
            }
        }

        for (bool s : supported) if (!s) return false;
        return true;
    }

    // The self-supporting printability constraints for per-vertex thickness
    // patterns take the form of inequality constraints on the position and
    // thickness parameters. These inequalities are linear apart from a min()
    // operation on the supporting candidates' z coordinates (all vertices are
    // supported if each is above the minimum of its neighbors.)
    // This min() is applied based on "params" so that a set of linear
    // inequality constraints is returned.
    // The constraints are in the form:
    //      C [p] >= 0
    //        [1]
    // where p is the parameter vector and C is the constraint matrix returned.
    //
    // For now, we assume all "dependency cycles" can be broken with a simple
    // heuristic from the pattern's default positions. By dependency cycles, we
    // mean we don't want both u to be a support candidate for v and v for u.
    //
    // The heuristic is as follows (**applied to default positions**)
    //   1) Mark all vertices supported from below with candidate lists.
    //   2) Mark vertices, v, unsupported from below with their neighbors that
    //      are supported by some other node than v.
    //      Assert there is only one for simplicity.
    //
    // For regular vertices, only allow a supporting vertex from below. For all
    // vertices that do not satisfy this constraint, allow a supporting vertex
    // from the same height, but assert that they are all supported vertices.
    Eigen::MatrixXd selfSupportingConstraints(
            const std::vector<Real> &params) const
    {
        if (thicknessType() != ThicknessType::Vertex)
            throw std::runtime_error("Only Per-Vertex Thickness Supported");

        std::vector<Edge> edges;
        std::vector<Point> points;
        std::vector<size_t> thicknessVars;
        std::vector<Eigen::Matrix3Xd> posMaps;
        printabilityGraph(params, points, edges, thicknessVars, posMaps);
        const size_t numPoints = points.size();

        // defaultZCoords: original z coordinate of each vertex in printability graph
        // currentZCoords: current    z coord of *bottom* of the sphere associated with each vtx
        // posMaps:        linear map expressing *bottom* of the sphere associated with each vtx
        std::vector<Real> currentZCoords, defaultZCoords;
        currentZCoords.reserve(numPoints), defaultZCoords.reserve(numPoints);
        for (size_t i = 0; i < numPoints; ++i) {
            defaultZCoords.push_back(m_printGraphVtx[i].pt[2]);
            // Determine sphere bottom by subtracting vertex radius from z coord
            currentZCoords.push_back(points[i][2] - params.at(thicknessVars[i]));
            posMaps[i](2, thicknessVars[i]) -= 1.0;
        }

        // Determine build platform height
        Real minZ = *std::min_element(defaultZCoords.begin(), defaultZCoords.end());

        // Build adjacency list
        std::vector<std::vector<size_t>> adj(numPoints);
        for (const auto &e : m_printGraphEdge) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }

        // Brute-force collection of the candidates to support a vertex.
        struct SupportCandidates {
            bool hasCandidate(size_t u) const {
                return std::find(candidates.begin(),
                                 candidates.end(), u) != candidates.end();
            }

            void add(size_t u) {
                if (!hasCandidate(u))
                    candidates.push_back(u);
            }

            void remove(size_t u) {
                auto it = std::find(candidates.begin(),
                                    candidates.end(), u);
                if (it == candidates.end()) throw std::runtime_error("Attempted to remove nonexistant support candidate: " + std::to_string(u));
                candidates.erase(it);
            }

            bool   empty() const { return candidates.empty(); }
            size_t count() const { return candidates.size(); }
            std::vector<size_t> candidates;
        };

        // Determine the candidates for supporting each vertex
        const double tol = m_tolerance();
        std::vector<SupportCandidates> supportCandidates(numPoints);
        std::vector<bool> needsSupport(numPoints, true);
        for (size_t u = 0; u < numPoints; ++u) {
            // Mark vertices on the build platform as not needing support
            if (defaultZCoords[u] < minZ + tol) {
                needsSupport[u] = false; continue;
            }
            for (size_t v : adj[u]) {
                // Must be definitively above, not horizontal.
                if (defaultZCoords[u] > defaultZCoords[v] + tol)
                    supportCandidates[u].add(v);
            }
        }

        // For each vertex without a supporting candidate, determine one of its
        // neighbors above that can be converted into a supporting candidate.
        for (size_t u = 0; u < numPoints; ++u) {
            if (!needsSupport[u]) continue;
            if (supportCandidates[u].empty()) {
                std::vector<size_t> options;
                for (size_t v : adj[u]) {
                    if (supportCandidates[v].hasCandidate(u)) {
                        // If there are other vertices than u that can support v,
                        // we can choose to make v support u.
                        if (supportCandidates[v].count() > 1) options.push_back(v);
                    }
                    else options.push_back(v);
                }
                if (options.size() == 0) {
                    std::cerr << u << " neighbors:";
                    for (size_t v : adj[u]) std::cerr << "\t" << v;
                    std::cerr << std::endl;
                    throw std::runtime_error("No options remain to support " + std::to_string(u));
                }

                // Take all options at the lowest available height.
                // TODO: check when this works in general.
                Real minOptionHeight = safe_numeric_limits<Real>::max();
                for (size_t v : options)
                    minOptionHeight = std::min(minOptionHeight, defaultZCoords[v]);
                for (size_t v : options) {
                    if (defaultZCoords[v] < minOptionHeight + tol) {
                        supportCandidates[u].add(v);
                        if (supportCandidates[v].hasCandidate(u)) {
                            // std::cerr << "removing " << u << " from " << v << std::endl;
                            supportCandidates[v].remove(u);
                        }
                    }
                }
            }
        }

        // Assert that we've found a valid set of support relationships.
        for (size_t u = 0; u < numPoints; ++u)
            assert(!(needsSupport[u] && supportCandidates[u].empty()));

#if 0
        for (size_t u = 0; u < numPoints; ++u) {
            if (!needsSupport[u]) {
                std::cerr << u + 1 << " doesn't need support" << std::endl;
                continue;
            }
            else {
                std::cerr << u + 1 << " supported by:";
                for (size_t v : supportCandidates[u].candidates)
                    std::cerr << "  " << v + 1;
                std::cerr << std::endl;
            }
        }
#endif
        // Create a constraint for every vertex that needs support
        std::vector<Eigen::Matrix<Real, 1, Eigen::Dynamic>> constraints;
        for (size_t u = 0; u < numPoints; ++u) {
            if (!needsSupport[u]) continue;
            // The constraint is that we stay above the lowest neighbor. We
            // apply the minimum operation here by looking at currentZCoords.
            const auto &sc = supportCandidates[u];
            const size_t NONE = std::numeric_limits<size_t>::max();
            size_t lowestCandidate = NONE;
            Real minHeight = safe_numeric_limits<Real>::max();
            for (size_t c : sc.candidates) {
                if (currentZCoords[c] < minHeight) {
                    minHeight = currentZCoords[c];
                    lowestCandidate = c;
                }
            }
            assert(lowestCandidate != NONE);
            // Constraint of the form c_i [p] >= 0
            //                            [1]
            constraints.push_back(posMaps[u].row(2) -
                                  posMaps[lowestCandidate].row(2));
            if (size_t(constraints.back().cols()) != params.size() + 1) {
                std::cerr << constraints.back().cols() << " cols and "
                          << params.size() << " parameters" << std::endl;
                std::cerr << posMaps[u].row(2) << std::endl;
                std::cerr << posMaps[lowestCandidate].row(2) << std::endl;
            }
            assert(size_t(constraints.back().cols()) == params.size() + 1);
        }

        // Detect and remove "constant" constraints not acting on variables;
        // these must be satisfied
        std::vector<Eigen::Matrix<Real, 1, Eigen::Dynamic>> prunedConstraints;
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
                std::cerr << "Detected constant constraint" << std::endl;
                if (c[params.size()] < -tol)
                    throw std::runtime_error("Infeasible: constant constraint unsatisfied");
            }
        }

        // TODO: think about reducing the constraint system

        // Convert constraints into Eigen format.
        Eigen::MatrixXd C(prunedConstraints.size(), params.size() + 1);
        for (size_t i = 0; i < prunedConstraints.size(); ++i)
            C.row(i) = prunedConstraints[i];
        return C;
    }

    // for a given vertex index, return parameters related to it
    std::vector<int> pointToParametersIndices(Point point) {
        std::vector<int> result;
        int baseIdx = m_findBaseVertex(point);
        int indepIdx = m_indepVtxForBaseVtx[baseIdx];
        int dofs = m_baseVertexPositioners[indepIdx].numDoFs();

        for (unsigned i=0; i<dofs; i++)
            result.push_back(m_baseVertexVarOffsets[baseIdx].position+i);

        result.push_back(m_baseVertexVarOffsets[baseIdx].thickness);
        result.push_back(m_baseVertexVarOffsets[baseIdx].blending);

        return result;
    }

    // for a given vertex index, return parameters related to it
    Point parameterIndexToPoint(size_t p) {
        // pass to original parameters index
        size_t correspondingOriginalIndex = m_filteredParamsToOriginal[p];

        Point result;

        if (correspondingOriginalIndex < numPositionParams()) {
            for (size_t i = 0; i < m_baseVertices.size(); ++i) {
                size_t numDofs = m_baseVertexPositioners[i].numDoFs();
                size_t offset = m_baseVertexVarOffsets[i].position;
                if (correspondingOriginalIndex >= offset && correspondingOriginalIndex < (offset + numDofs)) {
                    result = m_baseVertices[i];
                }
            }
        }
        else {
            for (size_t i = 0; i < m_baseVertices.size(); ++i) {
                if (correspondingOriginalIndex == m_baseVertexVarOffsets[i].thickness || correspondingOriginalIndex == m_baseVertexVarOffsets[i].blending) {
                    result = m_baseVertices[i];
                }
            }
        }

        return result;
    }


private:
    std::vector<decltype(PatternSymmetry::nodePositioner(Point()))> m_baseVertexPositioners;

    // Get the positioned graph and decoded per-vertex variables for a
    // particular graph (e.g. inflation graph) under parameter values "params"
    template<typename Real>
    void m_getGraphForParameters(
                        const std::vector<TransformedVertex> &graphVertices,
                        const std::vector<Edge>              &graphEdges,
                        const std::vector<Real> &params,
                        std::vector<Point3<Real>> &points,
                        std::vector<Edge> &edges,
                        std::vector<Real> &thicknesses,
                        std::vector<Real> &blendingParams) const {

        // Transform input parameters into original size
        std::vector<Real> inflationParams(m_originalParameters.size());
        std::vector<Real> originalParams(m_originalParameters.size());

        if (isAutodiffType<Real>())
            for (unsigned i = 0; i < originalParams.size(); i++)
                originalParams[i] = params[0]/params[0] * m_originalParameters[i]; //TODO: solve this without a trick
        else
            for (unsigned i = 0; i < originalParams.size(); i++)
                originalParams[i] = m_originalParameters[i];

        if (m_useFilterMask) {
            unsigned originalIdx=0;
            unsigned inputIdx=0;
            for (; originalIdx < originalParams.size(); originalIdx++) {
                inflationParams[originalIdx] = m_parametersMask[originalIdx] ? originalParams[originalIdx] : params[inputIdx++];
            }
        }
        else {
            inflationParams = params;
        }

        if (inflationParams.size() != numParams()) {
            std::cout << "Params size: " << inflationParams.size() << std::endl;
            std::cout << "Expected size: " << numParams() << std::endl;
            throw std::runtime_error("Invalid number of params.");
        }

        edges = graphEdges;

        // Position the base graph vertices using params
        std::vector<Point3<Real>> baseGraphPos;
        baseGraphPos.reserve(m_baseVertices.size());
        for (size_t i = 0; i < m_baseVertices.size(); ++i) {
            const auto &pos = m_baseVertexPositioners[i];
            size_t offset = m_baseVertexVarOffsets[i].position;
            baseGraphPos.push_back(pos.template getPosition<Real>(&inflationParams[offset]));
        }

        points.clear(), points.reserve(graphVertices.size());
        // Determine the vertex positions from the base positions.
        for (const TransformedVertex &iv : graphVertices)
            points.push_back(iv.iso.apply(baseGraphPos.at(iv.origVertex)));

        // Determine thicknesses and blending parameters
        if (thicknessType() == ThicknessType::Edge)
            throw std::runtime_error("Edge thickness currently unimplemented"); // will be easy to implement with m_inflEdgeOrigin

        thicknesses   .clear(), thicknesses   .reserve(graphVertices.size());
        blendingParams.clear(), blendingParams.reserve(graphVertices.size());

        for (const TransformedVertex &iv : graphVertices) {
            const auto &vo = m_baseVertexVarOffsets.at(iv.origVertex);
            thicknesses   .push_back(inflationParams.at(vo.thickness));
            blendingParams.push_back(inflationParams.at(vo.blending));
        }
    }

    virtual double m_tolerance() const override { return PatternSymmetry::tolerance; }

    // Find the base vertex within symmetry tolerance of p
    // (or throw an exception if none exists).
    size_t m_findBaseVertex(const Point &p) const;

    // How many edges to traverse outward when adding meshing-cell-adjacent
    // vertices to the inflation graph.
    size_t m_inflationNeighborhoodEdgeDist;
};

#include "WireMesh.inl"

#endif /* end of include guard: WIREMESH_HH */
