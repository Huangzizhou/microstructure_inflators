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
    for (int i = 0; i < numQuads; ++numQuads) {
        auto entry = params[i];
        m_allTopologies[i] = load_wire_mesh(entry["symmetry"], entry["pattern"]);
        m_allParameters[i] = entry["params"].get<std::vector<double>>();
    }

    m_V.resize(V.size(), 3);
    m_F.resize(F.size(), 4);
}

////////////////////////////////////////////////////////////////////////////////

// Set active quad in the background mesh
void WireQuadMesh::setActiveQuad(int idx) {
    //
}

// Retrieve inflation parameters for the active quad
std::vector<double> WireQuadMesh::params() const {
    std::vector<double> params;
    return params;
}

// Retrieve Jacobian matrix mapping the reference square [-1,1]² to the active quad
Eigen::Matrix3d WireQuadMesh::jacobian() const {
    return Eigen::Matrix3d();
}

// -----------------------------------------------------------------------------

void WireQuadMesh::inflationGraph(const std::vector<double> &allParams,
    std::vector<WireQuadMesh::Point> &stitchedPoints,
    std::vector<WireQuadMesh::Edge> &stitchedEdges,
    std::vector<double> &stitchedThicknesses,
    std::vector<double> &stitchedBlendingParams) const
{

}
