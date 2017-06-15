#include "InflatorTypes.hh"
#include "CGALClippedVolumeMesher.hh"
// #include "IGLSurfaceMesherMC.hh"
#include "PaperVisualizationSDFunc.hh"
#include "SignedDistance.hh"

int main(int argc, const char *argv[]) {
    if (argc != 7) {
        std::cerr << "usage: ./VisualizeQuantities r1 r2 r3 r4 s facet_distance" << std::endl;
        std::cerr << "example: ./VisualizeQuantities 0.25 0.25 0.25 0.25 0.04 2e-4" << std::endl;
        exit(-1);
    }

    double r1 = std::stod(argv[1]);
    double r2 = std::stod(argv[2]);
    double r3 = std::stod(argv[3]);
    double r4 = std::stod(argv[4]);
    double s  = std::stod(argv[5]);
    double facet_distance = std::stod(argv[6]);

    PaperVisualizationSDFunc sdfunc(r1, r2, r3, r4, s);

    // IGLSurfaceMesherMC mesher;
    CGALClippedVolumeMesher mesher;

    mesher.meshingOptions.facetDistance = facet_distance;
    mesher.meshingOptions.cellSize = 1.0;

    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;

    {
        std::vector<MeshIO::IOVertex > unionVertices;
        std::vector<MeshIO::IOElement> unionElements;

        sdfunc.mode = PaperVisualizationSDFunc::Mode::SHORT_EDGE1;
        mesher.mesh(sdfunc, unionVertices,  unionElements);

        sdfunc.mode = PaperVisualizationSDFunc::Mode::EDGE2;
        mesher.mesh(sdfunc, vertices,  elements);
        size_t offset = unionVertices.size();
        for (const auto &v : vertices) unionVertices.push_back(v);
        for (const auto &e : elements) unionElements.emplace_back(e[0] + offset, e[1] + offset, e[2] + offset, e[3] + offset);

        sdfunc.mode = PaperVisualizationSDFunc::Mode::EDGE3;
        mesher.mesh(sdfunc, vertices,  elements);
        offset = unionVertices.size();
        for (const auto &v : vertices) unionVertices.push_back(v);
        for (const auto &e : elements) unionElements.emplace_back(e[0] + offset, e[1] + offset, e[2] + offset, e[3] + offset);
        MeshIO::save("hard_union_joint1.msh", unionVertices, unionElements);
    }

    sdfunc.mode = PaperVisualizationSDFunc::Mode::BLEND_FULL;
    mesher.mesh(sdfunc, vertices,  elements);
    MeshIO::save("full_blend.msh", vertices, elements);

    sdfunc.mode = PaperVisualizationSDFunc::Mode::BLEND_HULL;
    mesher.mesh(sdfunc, vertices,  elements);
    MeshIO::save("hull_blend.msh", vertices, elements);

    sdfunc.mode = PaperVisualizationSDFunc::Mode::HULL;
    mesher.mesh(sdfunc, vertices,  elements);
    MeshIO::save("blending_region.msh", vertices, elements);

    {
        std::vector<MeshIO::IOVertex > unionVertices;
        std::vector<MeshIO::IOElement> unionElements;

        sdfunc.mode = PaperVisualizationSDFunc::Mode::JOINT_1;
        mesher.mesh(sdfunc, unionVertices,  unionElements);
        sdfunc.mode = PaperVisualizationSDFunc::Mode::JOINT_2;
        mesher.mesh(sdfunc, vertices,  elements);
        size_t offset = unionVertices.size();
        for (const auto &v : vertices) unionVertices.push_back(v);
        for (const auto &e : elements) unionElements.emplace_back(e[0] + offset, e[1] + offset, e[2] + offset, e[3] + offset);
        MeshIO::save("joint_pair_union.msh", unionVertices, unionElements);
    }

    sdfunc.mode = PaperVisualizationSDFunc::Mode::JOINT_UNION;
    mesher.mesh(sdfunc, vertices,  elements);
    MeshIO::save("joint_pair_smooth_union.msh", vertices, elements);

    sdfunc.mode = PaperVisualizationSDFunc::Mode::REDUCED_SMOOTH;
    mesher.mesh(sdfunc, vertices,  elements);
    MeshIO::save("reduced_smooth_union.msh", vertices, elements);

    sdfunc.mode = PaperVisualizationSDFunc::Mode::REDUCED_SMOOTH_BULGE;
    mesher.mesh(sdfunc, vertices,  elements);
    MeshIO::save("reduced_smooth_bulge.msh", vertices, elements);

    return 0;
}