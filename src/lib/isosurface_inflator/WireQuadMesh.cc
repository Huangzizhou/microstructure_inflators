#if HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#include "WireQuadMesh.hh"
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

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

WireQuadMesh::WireQuadMesh(
    const std::vector<MeshIO::IOVertex> &V,
    const std::vector<MeshIO::IOElement> &F,
    const nlohmann::json &params)
{
    size_t numQuads = F.size();
    m_allTopologies.clear(); m_allTopologies.resize(numQuads);
    m_allParameters.clear(); m_allParameters.resize(numQuads);
    m_allJacobians.clear(); m_allJacobians.resize(numQuads);
    for (size_t i = 0; i < numQuads; ++i) {
        auto entry = params[i];
        m_allParameters[i] = entry["params"].get<std::vector<double>>();
        m_allTopologies[i] = load_wire_mesh(entry["symmetry"], entry["pattern"]);
        assert(m_allTopologies[i]->thicknessType() == m_thicknessType);

        // Read jacobian as a row-major matrix, but Eigen matrices default
        // to column-major storage, hence the transpose
        // std::vector<double> x = entry["jacobian"];
        // Eigen::Matrix3d jacobian;
        // std::copy_n(x.data(), 9, jacobian.data());
        // jacobian.transposeInPlace();
        // m_allJacobians[i] = jacobian.transpose();
    }

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
}

////////////////////////////////////////////////////////////////////////////////

// Inflation parameters for the whole mesh (simple concatenation)
std::vector<double> WireQuadMesh::params() const {
    std::vector<double> params;
    for (const auto &p : m_allParameters) {
        params.insert(params.end(), p.begin(), p.end());
    }
    return params;
}

// Representative cell bounding box (region to be meshed)
BBox<Point3D> WireQuadMesh::boundingBox() const {
    Eigen::RowVector3d minV = m_V.colwise().minCoeff().array();
    Eigen::RowVector3d maxV = m_V.colwise().maxCoeff().array();
    Point3D a = minV.transpose();
    Point3D b = maxV.transpose();
    BBox<Point3D> bbox(a, b);
    return bbox;
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
    for (const auto &wmesh : m_allTopologies) {
        const size_t np = wmesh->numParams();
        std::vector<Real> wmparams(&allParams[paramOffset], &allParams[paramOffset] + np);
        paramOffset += np;
        std::vector<Point> points;
        std::vector<Edge> edges;
        std::vector<double> thicknesses;
        std::vector<double> blendingParams;
        wmesh->periodCellGraph(wmparams, points, edges, thicknesses, blendingParams);

        for (const auto &v : points)    allVerts.emplace_back(v);
        for (double t : thicknesses)    allThicknesses.push_back(t);
        for (double b : blendingParams) allBlendingParams.push_back(b);
    };

    // The vertices that were merged into a particular stitched output vertex.
    std::vector<std::vector<size_t>> m_stitchedVertices;
    std::vector<Edge> m_stitchedEdges;

    /////////////////////////////////////
    // TODO: Stitch graph on quad mesh //
    /////////////////////////////////////

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

    // _OutputGraph("test_inflation_graph.obj", stitchedPoints, stitchedEdges);
}

#endif
