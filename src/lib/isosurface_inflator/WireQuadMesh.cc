////////////////////////////////////////////////////////////////////////////////
#include "WireQuadMesh.hh"
////////////////////////////////////////////////////////////////////////////////

WireQuadMesh::WireQuadMesh(
    const std::vector<MeshIO::IOVertex> &V,
    const std::vector<MeshIO::IOElement> &F,
    const nlohmann::json &params)
{
    // Stuff
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
