#include "StitchedWireMesh.hh"
#include "PatternSignedDistance.hh"
#include "MidplaneMesher.hh"
#include "IGLSurfaceMesherMC.hh"

#include <memory>
#include <iostream>
#include <vector>
#include <string>

#include <cstdlib>
#include <ctime>

using namespace std;

using WMeshSPtr = std::shared_ptr<WireMeshBase>;

template<size_t N> struct DebugSymmetry    { using type = Symmetry::Square<>; };
template<>         struct DebugSymmetry<3> { using type = Symmetry::Cubic<>; };
template<size_t N> using DebugSymmetry_t = typename DebugSymmetry<N>::type;

template<size_t N>
void execute(const std::vector<std::string> &topologyPaths) {
    // Load all topologies passed on command line
    std::vector<WMeshSPtr> topologies;
    for (const auto &path : topologyPaths)
        topologies.emplace_back(make_shared<WireMesh<DebugSymmetry_t<N>>>(path));

    // Assign random topologies and parameters to each cell
    NDCubeArray<WMeshSPtr, N, 3> topologyGrid;
    NDCubeArray<std::vector<double>, N, 3> parameterGrid;

    topologyGrid.visit([&](WMeshSPtr &wm, const NDArrayIndex<N> &idxs) {
        wm = topologies[rand() % topologies.size()];
        parameterGrid(idxs) = wm->defaultParameters();
    });

    auto swm = make_stitched_wire_mesh(topologyGrid);

    auto params = swm.paramsFromParamGrid(parameterGrid);

    PatternSignedDistance<double, StitchedWireMesh<N>> sdf(swm);

    // Note: JointBlendMode could be set differently in MeshingOptions
    sdf.setParameters(params, Eigen::Matrix3d(), JointBlendMode::HULL);

    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;

    if (N == 2) MidplaneMesher().mesh(sdf, vertices, elements);
    if (N == 3) IGLSurfaceMesherMC().mesh(sdf, vertices, elements);
    MeshIO::save("meshed_cell.msh", vertices, elements);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "usage: TestStitchedWireMesh topology1.wire [topology2.wire...]" << std::endl;
        exit(-1);
    }

    srand(time(NULL));

    std::vector<std::string> topologyPaths;
    for (int i = 1; i < argc; ++i)
        topologyPaths.emplace_back(argv[i]);

    // Deduce dimension from the first topology
    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;
    MeshIO::load(topologyPaths[0], vertices, elements);
    BBox<Point3D>bb(vertices);

    if (bb.dimensions()[2] < 1e-5) execute<2>(topologyPaths);
    else                           execute<3>(topologyPaths);

    return 0;
}
