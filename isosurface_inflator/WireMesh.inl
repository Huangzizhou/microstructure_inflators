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
    m_adjacentVertices.clear(), m_adjacentEdges.clear();
    m_adjacentEdgeOrigin.clear();

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

    // Determine tiled vertices/edges adjacent to the base unit subgraph:
    // All non-identity isometries in the symmetry group can map some edges
    // outside the base unit while keeping them incident on the base unit. In
    // other words, one endpoint is mapped outside, while the other is mapped to
    // coincide with a base unit vertex).
    auto symmetryGroup = PatternSymmetry::symmetryGroup();
    std::vector<Point> adjacentMappedPoints;
    for (size_t ei = 0; ei < m_baseEdges.size(); ++ei) {
        const auto &e = m_baseEdges[ei];

        // TODO: add constraints for the periodic case
        size_t u = e.first, v = e.second;
        auto pu = m_baseVertices.at(u), pv = m_baseVertices.at(v);
        for (const auto &isometry : symmetryGroup) {
            if (isometry.isIdentity()) continue;
            auto mappedPu = isometry.apply(pu), mappedPv = isometry.apply(pv);
            bool puMappedInside = PatternSymmetry::inBaseUnit(mappedPu),
                 pvMappedInside = PatternSymmetry::inBaseUnit(mappedPv);
            Point outsideBasePoint;
            // Only care about edges where exactly one endpoint remains inside.
            if (puMappedInside != pvMappedInside) {
                size_t insideBaseIndex, outsideBaseIndex;
                if (puMappedInside) {
                    // Interior nodes should never stay inside!
                    assert(PatternSymmetry::nodeType(pu) != Symmetry::NodeType::Interior);
                    insideBaseIndex = m_findBaseVertex(mappedPu);
                    assert(insideBaseIndex == u); // will fail in TriplyPeriodic case
                    outsideBaseIndex = v;
                    outsideBasePoint = mappedPv;
                }
                if (pvMappedInside) {
                    // Interior nodes should never stay inside!
                    assert(PatternSymmetry::nodeType(pv) != Symmetry::NodeType::Interior);
                    insideBaseIndex = m_findBaseVertex(mappedPv);
                    assert(insideBaseIndex == v); // will fail in TriplyPeriodic case
                    outsideBaseIndex = u;
                    outsideBasePoint = mappedPu;
                }

                // Vertices and edges may be mapped to the same adjacent
                // vertex/edge multiple times by different isometries. For
                // example, if the z coordinate is zero, z reflection has no
                // effect. We detect these duplicates before they are created...
                // Determine if the point mapped outside is a duplicate.
                size_t adjacentPointIndex;
                auto api = std::find_if(adjacentMappedPoints.begin(), adjacentMappedPoints.end(),
                        [&](const Point &p) { return (p - outsideBasePoint).norm() < PatternSymmetry::tolerance; });
                if (api != adjacentMappedPoints.end()) {
                    adjacentPointIndex = std::distance(adjacentMappedPoints.begin(), api);
                    // std::cout << "Base points " << outsideBaseIndex  << ", "
                    //           << m_adjacentVertices[adjacentPointIndex].first
                    //           << " mapped to coincide under isometries "
                    //           << isometry << ", "
                    //           << m_adjacentVertices[adjacentPointIndex].second
                    //           << std::endl;
                }
                else {
                    adjacentPointIndex = m_adjacentVertices.size();
                    m_adjacentVertices.push_back({outsideBaseIndex, isometry});
                    adjacentMappedPoints.push_back(outsideBasePoint);
                }
                // Determine if the edge is a duplicate
                Edge e{insideBaseIndex, adjacentPointIndex};
                auto it = std::find(m_adjacentEdges.begin(), m_adjacentEdges.end(), e);
                if (it == m_adjacentEdges.end()) {
                    m_adjacentEdgeOrigin.push_back(ei);
                    m_adjacentEdges.push_back({insideBaseIndex, adjacentPointIndex});
                }
                else {
                    // size_t orig = std::distance(m_adjacentEdges.begin(), it);
                    // std::cout << "Edge " << e.first << ", "
                    //           << e.second << " is a duplicate "
                    //           << "(orig: " << m_adjacentEdgeOrigin[orig]
                    //           << ", current: " << ei << ")" << std::endl;
                }
            }
        }
    }

    // Enumerate position parameters:
    // Position parameters for each base vertex are determined using the
    // NodePositioner.
    m_baseVertexPositioners.reserve(m_baseVertices.size());
    for (const auto &p : m_baseVertices)
        m_baseVertexPositioners.push_back(PatternSymmetry::nodePositioner(p));
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

