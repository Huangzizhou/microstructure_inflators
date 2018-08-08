#if HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#include "WireQuadMesh.hh"
#include "MeshingOptions.hh"
#include "quadfoam/instantiate.h"
////////////////////////////////////////////////////////////////////////////////

namespace {

std::string lowercase(std::string data) {
    std::transform(data.begin(), data.end(), data.begin(), ::tolower);
    return data;
}

#define TRY_SYMMETRY(s, x, p)                                  \
    if (lowercase(x) == lowercase(#s))                         \
    {                                                          \
        return std::make_shared<WireMesh<Symmetry::s<>>>((p)); \
    }

#define TRY_KEY_VAL(s, a, x, p)                                \
    if (lowercase(x) == lowercase(#a))                         \
    {                                                          \
        return std::make_shared<WireMesh<Symmetry::s<>>>((p)); \
    }

std::shared_ptr<WireMeshBase> load_wire_mesh(const std::string &sym, const std::string &path) {
    TRY_SYMMETRY(Square, sym, path);
    TRY_SYMMETRY(Cubic, sym, path);
    TRY_SYMMETRY(Orthotropic, sym, path);
    TRY_SYMMETRY(Diagonal, sym, path);
    TRY_KEY_VAL(DoublyPeriodic, Doubly_Periodic, sym, path);
    TRY_KEY_VAL(TriplyPeriodic, Triply_Periodic, sym, path);
    return nullptr;
}

template<typename T>
void sort_unique(std::vector<T> &x) {
    std::sort(x.begin(), x.end());
    auto it = std::unique(x.begin(), x.end());
    x.resize(std::distance(x.begin(), it));
}


} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

WireQuadMesh::WireQuadMesh(
    const std::vector<MeshIO::IOVertex> &V,
    const std::vector<MeshIO::IOElement> &F,
    const nlohmann::json &params)
{
    // Read input parameters
    size_t numQuads = F.size();
    m_allTopologies.clear(); m_allTopologies.resize(numQuads);
    m_allParameters.clear(); m_allParameters.resize(numQuads);
    m_allJacobians.clear(); m_allJacobians.resize(numQuads);
    for (size_t i = 0; i < numQuads; ++i) {
        auto entry = params[i];
        m_allParameters[i] = entry["params"].get<std::vector<double>>();
        m_allTopologies[i] = load_wire_mesh(entry["symmetry"], entry["pattern"]);
        assert(m_allTopologies[i]->thicknessType() == m_thicknessType);
        m_allJacobians[i] = MeshingOptions::read_jacobian(entry["jacobian"]);
    }

    // Copy background quad mesh
    m_V.resize(V.size(), 3);
    m_F.resize(F.size(), 4);
    for (int v = 0; v < (int) V.size(); ++v) {
        for (int c = 0; c < 3; ++c) {
            m_V(v, c) = V[v][c];
        }
    }
    for (int f = 0; f < (int) F.size(); ++f) {
        assert(F[f].size() == 4);
        for (int lv = 0; lv < 4; ++lv) {
            m_F(f, lv) = F[f][lv];
        }
    }

    // Compute vertex id offset in concatenated graph
    m_vertexOffset.resize(numQuads + 1);
    m_vertexOffset.setZero();
    std::vector<Eigen::MatrixXd> allVertices;
    std::vector<Edge> allEdges;
    for (size_t f = 0; f < numQuads; ++f) {
        std::vector<Point> verts;
        std::vector<Edge>  edges;
        m_allTopologies[f]->periodCellGraph(verts, edges);
        m_vertexOffset[f+1] = m_vertexOffset[f] + verts.size();
        Eigen::MatrixXd V(verts.size(), 3);
        for (int i = 0; i < verts.size(); ++i) {
            V.row(i) = verts[i].transpose();
        }
        // Remap from [-1,1] to [0,1]
        V = V.array() * 0.5 + 0.5;
        allVertices.push_back(V);
        for (const auto &e : edges) {
            allEdges.emplace_back(e.first + m_vertexOffset[f], e.second + m_vertexOffset[f]);
        }
    };

    // Stitch input graph and remove duplicate vertices
    Eigen::MatrixXi tmpF;
    bool success = quadfoam::instanciate_pattern(
        m_V, m_F,
        allVertices, {},
        m_graphReducedVertices,
        tmpF,
        true,
        PatternSymmetry::tolerance,
        &m_graphFullToReduced
    );
    if (!success) {
        throw std::runtime_error("Could not stitch the wiremesh together!");
    }

    // Build inverse map `m_stitchedVertices`, as well as the set `m_stitchedEdges`
    m_stitchedVertices.resize(m_graphReducedVertices.size(), {});
    for (size_t i = 0; i < m_graphFullToReduced.size(); ++i) {
        int j = m_graphFullToReduced[i];
        assert(j < m_graphReducedVertices.size());
        m_stitchedVertices[j].push_back(i);
    }
    for (auto &e : allEdges) {
        m_stitchedEdges.emplace_back(m_graphFullToReduced[e.first], m_graphFullToReduced[e.second]);
    }
    sort_unique(m_stitchedEdges);
}

////////////////////////////////////////////////////////////////////////////////
// Bilinear map for a quad (a,b,c,d):
// (u,v) --> a + u*v*(a - b + c - d) + u*(-a + b) + v*(-a + d)


// Set currently active quad
void WireQuadMesh::setActiveQuad(int idx) {
    m_activeQuad = idx;
}

// -----------------------------------------------------------------------------

// Build the inflation graph for the active quad mesh, stitching together adjacent
// nodes (averaging stitched points' locations, thicknesses, and blending params).
//
// @param[in]  allParams               { Inflation parameters for the whole mesh }
// @param[out] stitchedPoints          { Graph vertices positions }
// @param[out] stitchedEdges           { Graph edges indices }
// @param[out] stitchedThicknesses     { Graph edge thicknesses }
// @param[out] stitchedBlendingParams  { Graph vertices blending parameters }
//
void WireQuadMesh::inflationGraph(const std::vector<double> &allParams,
    std::vector<WireQuadMesh::Point> &stitchedPoints,
    std::vector<WireQuadMesh::Edge> &stitchedEdges,
    std::vector<double> &stitchedThicknesses,
    std::vector<double> &stitchedBlendingParams) const
{
    std::vector<Point>  allVerts;
    std::vector<Edge>   allEdges;
    std::vector<double> allThicknesses;
    std::vector<double> allBlendingParams;

    // Build period cell graph for each wire mesh separately
    int paramOffset = 0;
    for (int i = 0; i < m_F.rows(); ++i) {
        const auto wmesh = m_allTopologies[i];
        const size_t np = wmesh->numParams();
        std::vector<Real> wmparams(&allParams[paramOffset], &allParams[paramOffset] + np);
        paramOffset += np;
        std::vector<Point> points;
        std::vector<Edge> edges;
        std::vector<double> thicknesses;
        std::vector<double> blendingParams;
        wmesh->periodCellGraph(wmparams, points, edges, thicknesses, blendingParams);

        // Graph parameters (radii & blending) are expressed in the world space, i.e. after
        // mapping the reference square [-1,1]² to the target parallelogram.
        // In order to stitch adjacent pattern inscribed in different parallelograms,
        // we first need to scale the radii and blending params by sqrt(jac^-1(i))`,
        // where i is the index of the element (quad), and jac is the jacobian (linear mapping)
        // mapping the reference square [-1,1]² to the parallelogram of the element i.

        double scaling = std::sqrt(m_allJacobians[i].inverse().determinant());

        for (const auto &v : points)    allVerts.emplace_back(v);
        for (double t : thicknesses)    allThicknesses.push_back(scaling * t);
        for (double b : blendingParams) allBlendingParams.push_back(scaling * b);
    };

    // Extracted stitched graph from the replicated graph by averaging
    // merged vertices' values.
    stitchedEdges = m_stitchedEdges;
    stitchedPoints.clear();
    stitchedThicknesses.clear();
    stitchedBlendingParams.clear();
    const size_t nsv = m_stitchedVertices.size();
    stitchedPoints.reserve(nsv);
    stitchedThicknesses.reserve(nsv);
    stitchedBlendingParams.reserve(nsv);
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

        //////////////////////////////////
        // TODO: Size based on Jacobian //
        //////////////////////////////////

        stitchedPoints.push_back(pt / sv.size());
        stitchedThicknesses.push_back(t / sv.size());
        stitchedBlendingParams.push_back(b / sv.size());
    }

    // Now, for each vertices, scale back the (radius, blending) parameters based
    // on the jacobian of the bilinear map (that maps the ref square [-1,1]² to
    // the active quad).



    // _OutputGraph("test_inflation_graph.obj", stitchedPoints, stitchedEdges);
}

#endif
