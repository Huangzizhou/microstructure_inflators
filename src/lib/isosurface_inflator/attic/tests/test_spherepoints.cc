#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include "SpherePoints.hh"
#include "JointViewer/ConvexHull.hh"
#include "JointViewer/ConvexHullTriangulation.hh"

int main(int argc, const char *argv[]) {
    size_t numPts = 68000;
    if (argc == 2) { numPts = std::stoi(argv[1]); }

    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;

    generateSpherePoints(numPts, vertices);
    for (auto &v : vertices) { v.point *= 0.5; v.point[0] -= 2.0; }

    generateSpherePoints(numPts, vertices);
    for (auto &v : vertices) { v.point *= 0.5; v.point[1] -= 2.0; }

    generateSpherePoints(numPts, vertices);

    for (size_t i = 0; i < vertices.size() - 1; ++i) { elements.emplace_back(i, i + 1); }
    MeshIO::save("spherePts.msh", vertices, elements);

    std::vector<MeshIO::IOVertex > hullVertices;
    std::vector<MeshIO::IOElement> hullElements;

    BENCHMARK_START_TIMER("convex_hull_3 Version");
    convexHull(vertices, hullVertices, hullElements);
    BENCHMARK_STOP_TIMER("convex_hull_3 Version");

    BENCHMARK_START_TIMER("Triangulation Version");
    convexHullFromTriangulation(vertices, hullVertices, hullElements);
    BENCHMARK_STOP_TIMER("Triangulation Version");

    BENCHMARK_START_TIMER("Triangulation Version with index tracking");
    std::vector<size_t> origIndex;
    convexHullFromTriangulation(vertices, hullVertices, hullElements, origIndex);
    BENCHMARK_STOP_TIMER("Triangulation Version with index tracking");

    MSHFieldWriter writer("spherePtsHullIndexed.msh", hullVertices, hullElements);
    ScalarField<double> idx(origIndex.size());
    for (size_t i = 0; i < origIndex.size(); ++i) idx[i] = origIndex[i];
    writer.addField("idx", idx, DomainType::PER_NODE);

    assert(hullVertices.size() == origIndex.size());
    for (size_t i = 0; i < hullVertices.size(); ++i) {
        assert(hullVertices[i].point[0] == vertices[origIndex[i]].point[0]);
        assert(hullVertices[i].point[1] == vertices[origIndex[i]].point[1]);
        assert(hullVertices[i].point[2] == vertices[origIndex[i]].point[2]);
    }

    MeshIO::save("spherePtsHull.msh", hullVertices, hullElements);

    BENCHMARK_REPORT();

    return 0;
}
