////////////////////////////////////////////////////////////////////////////////
// StitchedWireMesh.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Construct the stitched inflation graph for a period cell and its
//  surrounding neighbors (to implement tilings with spatially varying
//  parameters). Stitched vertex positions and radii/blending params are
//  averaged.
//
//  All unspecified neighbors are replaced with a copy of the central cell (for
//  proper handling of the tiling borders).
//
//  Either 2N or (3^N - 1) neighbors can be included (i.e. only face-incident
//  cells or also including  corner- and edge-incident cells), as requested by
//  boolean flag "faceNeighborsOnly".
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Created:  07/17/2017 13:05:10
////////////////////////////////////////////////////////////////////////////////
#ifndef STITCHEDWIREMESH_HH
#define STITCHEDWIREMESH_HH

#include <Utilities/NDArray.hh>
#include <function_traits.hh>
#include <Geometry.hh>
#include <stdexcept>
#include <PeriodicBoundaryMatcher.hh> // for FaceMembership
#include <memory>
#include <set>
#include "Symmetry.hh"
#include "WireMesh.hh"

// Stencil Visitors:
//  Selectively visit cells of the grid if they belong to the neighborhood
//  stencil. Wraps two types of visitors:
//  1) single-arg functors (accessing value only), and
//  2) two-arg functors (taking value and an NDArrayIndex)
// N:   Dimension of grid to visit
// FNO: Visit Face-adjacent Neighbors Only (true) or all grid cells (false)
// V:   Visitor functor type (to be deduced)
template<size_t N, bool FNO, class V> struct StencilVisitor;

template<size_t N, bool FNO = true> class StitchedWireMesh;
template<size_t N, bool FNO = true> StitchedWireMesh<N, FNO>
make_stitched_wire_mesh(const NDCubeArray<std::shared_ptr<WireMeshBase>, N, 3> &wmeshGrid) {
    return StitchedWireMesh<N, FNO>(wmeshGrid);
}

template<size_t N, bool FNO>
class StitchedWireMesh {
public:
    using PatternSymmetry = Symmetry::Null<>;
    using WMeshSPtr = std::shared_ptr<WireMeshBase>;
    using PVec = std::vector<double>;

    using Point = WireMeshBase::Point;
    using Edge  = WireMeshBase::Edge;

    // Alias for making a "masked" stencil visitor from an ordinary value visitor.
    template<class V> static StencilVisitor<N, FNO, V> SV(V &&v) { return StencilVisitor<N, FNO, V>(std::forward<V>(v)); }

    StitchedWireMesh(const NDCubeArray<WMeshSPtr, N, 3> &wmeshGrid) {
        set(wmeshGrid);
    }

    // Build merged graph, tracking which global replicated vertices get merged
    // into each stitched output vertex
    void set(const NDCubeArray<WMeshSPtr, N, 3> &wmeshGrid) {
        m_wmeshGrid = wmeshGrid;

        // Fill in missing wire meshes with the center mesh.
        auto centerWMesh = m_wmeshGrid.getCenter();
        if (!centerWMesh) throw std::runtime_error("Center wire mesh must be present");
        m_wmeshGrid.visit(SV([=](WMeshSPtr &wm) { if (!wm) wm = centerWMesh; }));

        double cellSize = Symmetry::TriplyPeriodic<>::template representativeMeshCell<double>().dimensions()[0];
        // Get period cell meshes of all wire mesh topologies in the grid
        // (default position params); offset by the grid mesh size.
        std::vector<Point> allVerts;
        std::vector<Edge>  allEdges;
        m_wmeshGrid.visit(SV([&](const WMeshSPtr &wmesh, const NDArrayIndex<N> &idxs) {
            Point cellCornerOffset = m_getCellOffset(cellSize, idxs);
            size_t offset = allVerts.size();

            std::vector<Point> verts;
            std::vector<Edge>  edges;
            wmesh->periodCellGraph(verts, edges);
            for (const auto &e : edges) allEdges.push_back({e.first + offset, e.second + offset});
            for (const auto &v : verts) allVerts.emplace_back(v + cellCornerOffset);
        }));

        // Make sure all wire meshes are of the same thickness type.
        m_thicknessType = centerWMesh->thicknessType();
        m_wmeshGrid.visit(SV([=](const WMeshSPtr &wm) {
            if (wm->thicknessType() != m_thicknessType)
                throw std::runtime_error("Mismatched thickness type");
        }));

        if (m_thicknessType != ThicknessType::Vertex)
            throw std::runtime_error("Only per-vertex thickness currently supported");

        // Brute-force merge of vertices within PatternSymmetry::TOL of each other.
        // Tie each output vertex back to the set of vertices merging into it.
        std::vector<int> stitchedIndex(allVerts.size(), -1);
        m_stitchedVertices.clear();
        const double tolerance = PatternSymmetry::tolerance;
        for (size_t i = 0; i < allVerts.size(); ++i) {
            const auto &vi = allVerts[i];
            if (stitchedIndex[i] >= 0) continue;
            stitchedIndex[i] = m_stitchedVertices.size();
            m_stitchedVertices.emplace_back(1, i);
            for (size_t j = i + 1; j < allVerts.size(); ++j) {
                if (stitchedIndex[j] >= 0) continue;
                const auto &vj = allVerts[j];
                if ((vi - vj).squaredNorm() < tolerance * tolerance) {
                    m_stitchedVertices.back().push_back(j);
                    stitchedIndex[j] = stitchedIndex[i];
                }
            }
        }

        // Re-link edges in allEdges to the stitched vertices
        // Also, remove duplicate edges (if topologies have edges on interface)
        std::set<UnorderedPair> uniqueEdges;
        for (auto &e : allEdges) {
            uniqueEdges.emplace(stitchedIndex[e.first],
                                stitchedIndex[e.second]);
        }
        m_stitchedEdges.clear(), m_stitchedEdges.reserve(uniqueEdges.size());
        for (auto &e : uniqueEdges) m_stitchedEdges.emplace_back(e[0], e[1]);

        std::vector<PointND<N>> stitchedVtxLoc;
        for (const auto &sv : m_stitchedVertices)
            stitchedVtxLoc.push_back(truncateFrom3D<PointND<N>>(allVerts.at(sv.at(0))));

        if (!FNO) {
            // Validate merge: check that all dangling vertices are on the boundary.
            // Note: only works on full (non-plus) stencil
            BBox<PointND<N>> bb(stitchedVtxLoc);
            std::vector<size_t> valence(stitchedVtxLoc.size());
            for (const auto &e : m_stitchedEdges) {
                ++valence.at(e.first);
                ++valence.at(e.second);
            }
            for (size_t i = 0; i < valence.size(); ++i) {
                assert(valence[i] > 0);
                using FM = PeriodicBoundaryMatcher::FaceMembership<N>;
                if ((valence[i] == 1) && !FM(stitchedVtxLoc[i], bb).onAnyFace())
                    throw std::runtime_error("Dangling interior vertex in stitched graph; topology mismatch?");
            }
        }

        // Output stitched graph for debugging
        // _OutputGraph("test_stitch_graph.msh", stitchedVtxLoc, m_stitchedEdges);

        m_numParams = 0;
        m_wmeshGrid.visit(SV([&](const WMeshSPtr &wm) { m_numParams += wm->numParams(); }));
    }

    size_t numParams() const { return m_numParams; }

    // Parameter vector for parameter grid: simply concatenate...
    // (copy instead of link)
    std::vector<double> paramsFromParamGrid(NDCubeArray<PVec, N, 3> &paramGrid) const {
        // Overwrite any missing parameter vectors with the center parameters.
        const PVec &centerParams = paramGrid.getCenter();
        paramGrid.visit(SV([&centerParams](PVec &pvec) {
            if (pvec.size() == 0) pvec = centerParams;
        }));

        // Build single parameter vector by concatenating in scanline order.
        std::vector<double> params;
        params.reserve(numParams());
        paramGrid.visit(SV([&params](const PVec &pvec) { for (double p : pvec) params.push_back(p); }));
        if (params.size() != numParams()) throw std::runtime_error("Incorrect parameter vector size");
        return params;
    }

    // Build the period cell graph for every cell, then stitch
    // together (averaging stitched points' locations, thicknesses, and
    // blending params).
    void inflationGraph(const std::vector<double> &allParams,
                        std::vector<Point>  &stitchedPoints,
                        std::vector<Edge>   &stitchedEdges,
                        std::vector<double> &stitchedThicknesses,
                        std::vector<double> &stitchedBlendingParams) const {
        size_t paramOffset = 0;
        double cellSize = Symmetry::TriplyPeriodic<>::template representativeMeshCell<double>().dimensions()[0];
        if (allParams.size() != numParams()) throw std::runtime_error("Incorrect parameter vector size");

        std::vector<Point>  allVerts;
        std::vector<Edge>   allEdges;
        std::vector<double> allThicknesses;
        std::vector<double> allBlendingParams;

        // Build period cell graph for each wire mesh separately
        m_wmeshGrid.visit(SV([&](const WMeshSPtr &wmesh, const NDArrayIndex<N> &idxs) {
            const size_t np = wmesh->numParams();
            std::vector<Real> wmparams(&allParams[paramOffset], &allParams[paramOffset] + np);
            paramOffset += np;
            std::vector<Point> points;
            std::vector<Edge> edges;
            std::vector<double> thicknesses;
            std::vector<double> blendingParams;
            wmesh->periodCellGraph(wmparams, points, edges, thicknesses, blendingParams);

            Point cellCornerOffset = m_getCellOffset(cellSize, idxs);
            for (const auto &v : points)    allVerts.emplace_back(v + cellCornerOffset);
            for (double t : thicknesses)    allThicknesses.push_back(t);
            for (double b : blendingParams) allBlendingParams.push_back(b);
        }));

        // Extracted stitched graph from the replicated graph by averaging
        // merged vertices' values.
        stitchedEdges = m_stitchedEdges;
        stitchedPoints.clear(), stitchedThicknesses.clear(), stitchedBlendingParams.clear();
        const size_t nsv = m_stitchedVertices.size();
        stitchedPoints.reserve(nsv), stitchedThicknesses.reserve(nsv), stitchedBlendingParams.reserve(nsv);
        for (const auto &sv : m_stitchedVertices) {
            // Take mean of location, thickness, and blending of all vertices
            // merged into this stitched output vertex.
            Point pt(Point::Zero());
            double t = 0, b = 0.0;
            for (size_t v : sv) {
                pt += allVerts.at(v);
                t  += allThicknesses.at(v);
                b  += allBlendingParams.at(v);
            }
            stitchedPoints.push_back(pt / sv.size());
            stitchedThicknesses.push_back(t / sv.size());
            stitchedBlendingParams.push_back(b / sv.size());
        }

        // _OutputGraph("test_inflation_graph.msh", stitchedPoints, stitchedEdges);
    }

    // The stitching process can violate the self-supporting constraint.
    // However, assuming each individual cell is printable, we can restore
    // printability by clamping each repositioned vertex to be at the same
    // level as the adjacent vertex it needs to support.
    template<typename Real>
    void printableInflationGraph(const std::vector<Real> &/* params */,
                        std::vector<Point3<Real>> &/* points */,
                        std::vector<Edge> &/* edges */,
                        std::vector<Real> &/* thicknesses */,
                        std::vector<Real> &/* blendingParams */) const {
        throw std::runtime_error("printableInflationGraph unimplemented");
    }

    ThicknessType thicknessType() const { return m_thicknessType; }

private:
    NDCubeArray<WMeshSPtr, N, 3> m_wmeshGrid;
    size_t m_numParams = 0;
    std::vector<Edge> m_stitchedEdges;
    // The vertices that were merged into a particular stitched output vertex.
    std::vector<std::vector<size_t>> m_stitchedVertices;
    ThicknessType m_thicknessType;

    Point m_getCellOffset(double cellSize, const NDArrayIndex<N> &idxs) const {
        Point offset(Point::Zero());
        for (size_t i = 0; i < N; ++i)
            offset[i] = cellSize * (double(idxs[i]) - 1);
        return offset;
    }
};

////////////////////////////////////////////////////////////////////////////////
// Implementation detail
////////////////////////////////////////////////////////////////////////////////
namespace detail {
    template<class V, typename ValType, size_t N>
    auto apply(V &&visitor, ValType &&val, const NDArrayIndex<N> & /* idxs */) -> decltype(visitor(val)) {
        return visitor(std::forward<ValType>(val));
    }

    template<class V, typename ValType, size_t N>
    auto apply(V &&visitor, ValType &&val, const NDArrayIndex<N> &idxs) -> decltype(visitor(val, idxs)) {
        return visitor(std::forward<ValType>(val), idxs);
    }

    template<size_t N>
    int l1DistFromCenter(const NDArrayIndex<N> &idxs) {
        int dist = 0;
        for (size_t i = 0; i < N; ++i) { dist += std::abs(int(idxs[i]) - 1); }
        return dist;
    }
    template<size_t N>
    bool inPlusStencil(const NDArrayIndex<N> &idxs) { return l1DistFromCenter(idxs) < 2; }
}

template<size_t N, bool FNO, class V>
struct StencilVisitor {
    StencilVisitor(V &&v) : visitor(v) { }
    StencilVisitor(const StencilVisitor &b) : visitor(b.visitor) { }
    using ValType = typename function_traits<V>::template arg<0>::type;

    void operator()(ValType &&val, const NDArrayIndex<N> &idxs) {
        if (FNO && !detail::inPlusStencil(idxs)) return;
        detail::apply(visitor, std::forward<ValType>(val), idxs);
    }

    V &visitor;
};

#endif /* end of include guard: STITCHEDWIREMESH_HH */
