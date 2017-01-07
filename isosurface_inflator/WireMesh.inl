#include <algorithm>
#include <iterator>

// Set from embedded graph
template<ThicknessType thicknessType, class Sym>
void WireMesh<thicknessType, Sym>::
set(const std::vector<MeshIO::IOVertex > &inVertices,
    const std::vector<MeshIO::IOElement> &inElements)
{
    for (const auto &e : inElements)
        if (e.size() != 2) throw std::runtime_error("Expected line element mesh.");

    // Clear old state.
    m_fullVertices.clear(), m_fullEdges.clear();
    m_baseVertices.clear(), m_baseEdges.clear();
    m_baseVertexPositioners.clear();

    m_fullVertices.reserve(inVertices.size());
    m_fullEdges.reserve(inElements.size());
    for (const auto &e : inElements) {
        assert(e.size() == 2);
        m_fullEdges.push_back({e[0], e[1]});
    }

    // Convert vertices to Point type, scaling into [-1, 1]
    BBox<Point> bb(inVertices);
    auto dim = bb.dimensions();
    for (const auto &v : inVertices) {
        // transform graph to [-1, 1]
        Point p;
        p.setZero();
        for (size_t i = 0; i < 3; ++i) {
            // Hack to handle 2D case: leave "empty dimension" coordinates at zero.
            if (std::abs(dim[i]) < 1e-6) {
                assert(i == 2);
                continue;
            }
            p[i] = (v[i] - bb.minCorner[i]) * (2.0 / dim[i]) - 1.0;
        }
        m_fullVertices.push_back(p);
    }

    // Determine the vertices in the symmetry base unit subgraph
    std::vector<int> baseVertexIndex(m_fullVertices.size(), -1);
    for (size_t i = 0; i < m_fullVertices.size(); ++i) {
        if (PatternSymmetry::inBaseUnit(m_fullVertices[i])) {
            baseVertexIndex[i] = m_baseVertices.size();
            m_baseVertices.push_back(m_fullVertices[i]);
        }
    }

    // Compute edges in the induced subgraph
    for (const auto &e : m_fullEdges) {
        int u = baseVertexIndex.at(e.first), v = baseVertexIndex.at(e.second);
        if ((u >= 0) && (v >= 0))
            m_baseEdges.push_back({u, v});
    }

    // Enumerate position parameters:
    // Position parameters for each base vertex are determined using the
    // NodePositioner.
    m_baseVertexPositioners.reserve(m_baseVertices.size());
    for (const auto &p : m_baseVertices)
        m_baseVertexPositioners.push_back(PatternSymmetry::nodePositioner(p));

    ////////////////////////////////////////////////////////////////////////////
    // Construct inflation graph
    ////////////////////////////////////////////////////////////////////////////
    // Determine vertices/edges in the inflation graph
    std::vector<TransformedVertex> rVertices;
    std::vector<TransformedEdge>   rEdges;
    replicatedGraph(PatternSymmetry::symmetryGroup(), rVertices, rEdges);

    // Detect edges inside or incident base cell.
    std::vector<bool> touchingBaseCell;
    for (const auto &re : rEdges) {
        touchingBaseCell.push_back(PatternSymmetry::inBaseUnit(rVertices[re.e[0]].pt) ||
                                   PatternSymmetry::inBaseUnit(rVertices[re.e[1]].pt));
    }

    // Keep edges inside or incident base cell
    std::vector<bool> keepEdge = touchingBaseCell;

    if (KEEP_BC_ADJ_JOINTS) {
        // Also keep the edges for joints incident these edges.
        std::vector<bool> keepJoint(rVertices.size(), false);
        for (size_t ei = 0; ei < rEdges.size(); ++ei) {
            if (!keepEdge[ei]) continue;
            const auto &re = rEdges[ei];
            keepJoint.at(re.e[0]) = keepJoint.at(re.e[1]) = true;
        }

        for (size_t ei = 0; ei < rEdges.size(); ++ei) {
            const auto &re = rEdges[ei];
            if (keepJoint.at(re.e[0]) || keepJoint.at(re.e[1]))
                keepEdge[ei] = true;
        }
    }

    // Copy over kept edges.
    m_inflEdge.clear(), m_inflEdgeOrigin.clear();
    std::vector<bool> vtxSeen(rVertices.size(), false);
    for (size_t ei = 0; ei < rEdges.size(); ++ei) {
        const auto &re = rEdges[ei];
        if (!keepEdge[ei]) continue;
        size_t u = re.e[0],
               v = re.e[1];
        m_inflEdge.push_back({u, v});
        m_inflEdgeOrigin.push_back(re.origEdge);
        vtxSeen.at(u) = vtxSeen.at(v) = true;
    }

    // Only keep vertices belonging to kept edges.
    m_inflVtx.clear();
    std::vector<size_t> vertexRenumber(rVertices.size(), std::numeric_limits<size_t>::max());
    for (size_t i = 0; i < rVertices.size(); ++i) {
        if (!vtxSeen[i]) continue;
        vertexRenumber[i] = m_inflVtx.size();
        m_inflVtx.push_back(rVertices[i]);
    }

    // Reindex kept edges.
    for (auto &e : m_inflEdge) {
        e.first = vertexRenumber.at(e.first);
        e.second = vertexRenumber.at(e.second);
        assert((e.first < m_inflVtx.size()) && (e.second < m_inflVtx.size()));
    }

    ////////////////////////////////////////////////////////////////////////////
    // Construct printability graph
    ////////////////////////////////////////////////////////////////////////////
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

    std::vector<TransformedEdge> pEdges;
    replicatedGraph(printabiltyIsometries,
                    m_printGraphVtx, pEdges);
    m_printGraphEdge.clear(), m_printGraphEdge.reserve(pEdges.size());
    for (const auto &pe : pEdges)
        m_printGraphEdge.push_back({pe.e[0], pe.e[1]});
}

// Set from embedded graph stored in obj/msh format.
template<ThicknessType thicknessType, class Sym>
void WireMesh<thicknessType, Sym>::
load(const std::string &wirePath) {
    std::vector<MeshIO::IOVertex> inVertices;
    std::vector<MeshIO::IOElement> inElements;
    MeshIO::load(wirePath, inVertices, inElements);
    set(inVertices, inElements);
}

template<class Point>
void _OutputGraph(const std::string &path, const std::vector<Point> &points,
                  const std::vector<std::pair<size_t, size_t>> &edges) {
    std::vector<MeshIO::IOVertex>  outVertices;
    std::vector<MeshIO::IOElement> outElements;
    outVertices.reserve(points.size()), outElements.reserve(edges.size());
    for (const Point &p : points) outVertices.emplace_back(p);
    for (const auto  &e :  edges) outElements.emplace_back(e.first, e.second);
    MeshIO::save(path, outVertices, outElements);
}

// Save the full embedded graph in obj/msh format.
template<ThicknessType thicknessType, class Sym>
void WireMesh<thicknessType, Sym>::
save(const std::string &path) const {
    _OutputGraph(path, m_fullVertices, m_fullEdges);
}

// Save the embedded base unit graph in obj/msh format.
template<ThicknessType thicknessType, class Sym>
void WireMesh<thicknessType, Sym>::
saveBaseUnit(const std::string &path) const {
    _OutputGraph(path, m_baseVertices, m_baseEdges);
}

// For symmetry debugging:
// save the tiled base unit graph in obj/msh format
template<ThicknessType thicknessType, class Sym>
void WireMesh<thicknessType, Sym>::
saveReplicatedBaseUnit(const std::string &path) const {
    std::vector<MeshIO::IOVertex> outVertices;
    std::vector<MeshIO::IOElement> outEdges;
    for (const auto &iso : PatternSymmetry::symmetryGroup()) {
        size_t offset = outVertices.size();
        for (const Point &p : m_baseVertices) { outVertices.emplace_back(iso.apply(p)); }
        for (const Edge  &e :    m_baseEdges) { outEdges.emplace_back(e.first + offset, e.second + offset); }
    }
    MeshIO::save(path, outVertices, outEdges);
}

// Find the base vertex within symmetry tolerance of p
// (or throw an exception if none exists).
// The number of base vertices to check is generally small, so fancy
// data structures shouldn't be needed.
template<ThicknessType thicknessType, class Sym>
size_t WireMesh<thicknessType, Sym>::
m_findBaseVertex(const Point &p) const {
    for (size_t i = 0; i < m_baseVertices.size(); ++i) {
        if ((p - m_baseVertices[i]).norm() < PatternSymmetry::tolerance)
            return i;
    }
    std::cout << "Failed to find " << p[0] << " " << p[1] << " " << p[2] << std::endl;
    for (const Point &v : m_baseVertices)
        std::cout << "    candidate " << v[0] << " " << v[1] << " " << v[2] << std::endl;
    throw std::runtime_error("Couldn't find corresponding base vertex.");
}

// For symmetry debugging:
// save the inflation graph in obj/msh format
template<ThicknessType thicknessType, class Sym>
void WireMesh<thicknessType, Sym>::
saveInflationGraph(const std::string &path, std::vector<double> params) const {
    if (params.size() == 0)
        params = defaultParameters();
    std::vector<Edge> igraphEdges;
    std::vector<Point> igraphVertices;
    std::vector<double> thicknesses, blendingParams;
    inflationGraph(params, igraphVertices, igraphEdges, thicknesses, blendingParams);

    _OutputGraph(path, igraphVertices, igraphEdges);
}

