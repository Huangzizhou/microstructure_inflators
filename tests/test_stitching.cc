#if HAS_LIBIGL

////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/MidplaneMesher.hh>
#include <isosurface_inflator/IGLSurfaceMesherMC.hh>
#include <isosurface_inflator/WireQuadMesh.hh>
#include <MeshFEM/StringUtils.hh>
// #include <catch2/catch.hpp>
#include <memory>
////////////////////////////////////////////////////////////////////////////////

const nlohmann::json default_meshing_options = R"({
    "domainErrorBound"        : 1e-5,
    "facetAngle"              : 30.0,
    "facetSize"               : 0.02,
    "facetDistance"           : 1e-3,
    "cellSize"                : 0.50,
    "edgeSize"                : 0.015,
    "cellRadiusEdgeRatio"     : 2.0,

    "marchingSquaresGridSize" : 128,
    "marchingSquaresCoarsening": 1,
    "marchingCubesGridSize"   : 128,

    "forceConsistentInterfaceMesh": false,
    "forceMSGridSize": true,

    "maxArea"                 : 0.001,
    "featureAngleThreshold"   : 0.7853981633974483
})"_json;

using json = nlohmann::json;

////////////////////////////////////////////////////////////////////////////////

template<typename Symmetry> struct SymmetryTraits { };

#define SYMMETRY_NAME(sym, name)                  \
    template<>                                    \
    struct SymmetryTraits<sym> {                  \
        static constexpr char value[] = name;     \
    };                                            \
    constexpr char SymmetryTraits<sym>::value []; \

SYMMETRY_NAME(Symmetry::Cubic<>, "cubic")
SYMMETRY_NAME(Symmetry::Orthotropic<>, "orthotropic")
SYMMETRY_NAME(Symmetry::Diagonal<>, "diagonal")
SYMMETRY_NAME(Symmetry::Square<>, "square")
SYMMETRY_NAME(Symmetry::TriplyPeriodic<>, "triply_periodic")
SYMMETRY_NAME(Symmetry::DoublyPeriodic<>, "doubly_periodic")
SYMMETRY_NAME(Symmetry::NonPeriodic<>, "non_periodic")

////////////////////////////////////////////////////////////////////////////////

template<int N, typename SymmetryType>
void test_quad_mesh(const std::string &mesh, const std::string &topology) {
    // Create mesher and load meshing options
    std::unique_ptr<MesherBase> mesher;
    if (N == 2) { mesher = std::make_unique<MidplaneMesher>(); }
    if (N == 3) { mesher = std::make_unique<IGLSurfaceMesherMC>(); }
    mesher->meshingOptions.load(default_meshing_options);

    // Load input quad mesh
    std::vector<MeshIO::IOVertex> vertices_in;
    std::vector<MeshIO::IOElement> elements_in;
    auto mesh_type = MeshIO::load(mesh, vertices_in, elements_in);
    if (mesh_type != MeshIO::MESH_QUAD) {
        std::cerr << "Invalid element type, expected quads." << std::endl;
        return;
    }

    // Build dummy json input
    json data = json::array();
    {
        std::vector<MeshIO::IOVertex> VI;
        std::vector<MeshIO::IOElement> FI;
        MeshIO::load(topology, VI, FI);
        WireMesh<SymmetryType> wm(VI, FI);
        for (int i = 0; i < (int) elements_in.size(); ++i) {
            json entry;
            entry["params"] = wm.defaultParameters();
            entry["symmetry"] = SymmetryTraits<SymmetryType>::value;
            entry["pattern"] = topology;
            data.push_back(entry);
        }
    }

    // Wire mesh embedded into a quad mesh
    WireQuadMesh wm(vertices_in, elements_in, data);

    for (int index = -1; index < 0 /* (int) elements_in.size() */; ++index) {
        // Set SDF
        wm.setActiveQuad(index);
        PatternSignedDistance<double, WireQuadMesh, WireQuadMesh::MapToBaseUnit> sdf(wm);
        sdf.setParameters(wm.params(), mesher->meshingOptions.jacobian, mesher->meshingOptions.jointBlendingMode);
        sdf.setMapFunctor(wm.mapFunctor());
        sdf.setBoundingBox(wm.boundingBox());

        // Mesh pattern and save result
        {
            std::vector<MeshIO::IOVertex> vertices_out;
            std::vector<MeshIO::IOElement> elements_out;

            mesher->meshInterfaceConsistently = true;
            mesher->mesh(sdf, vertices_out, elements_out);

            #ifdef DUMP_OUTPUT
            std::string basename = MeshFEM::replace_ext(MeshFEM::split(mesh, "/").back(), "");
            std::string suffix = (index < 0 ? "whole" : "q" + std::to_string(index));
            MeshIO::save(basename + "_" + SymmetryTraits<SymmetryType>::value + "_" + suffix + ".obj", vertices_out, elements_out);
            #endif
        }
    }
}

// -----------------------------------------------------------------------------

template<typename SymmetryType>
void test_quad_meshes(const std::vector<std::string> &meshes, const std::string &pattern) {
    for (auto filename : meshes) {
        test_quad_mesh<2, SymmetryType>(filename, pattern);
    }
}

////////////////////////////////////////////////////////////////////////////////

#define SECTION(x)

// TEST_CASE("inflate_and_stitch", "[isosurface_inflation]") {
int main(void) {
    std::string pattern_2d = DATA_DIR "patterns/2D/topologies/0001.obj";

    std::vector<std::string> meshes = {
        DATA_DIR "tests/quad_grid_orient_fuzzy.obj",
        DATA_DIR "tests/quad_grid_orient_perfect.obj",
        DATA_DIR "tests/quad_irregular_orient_fuzzy.obj",
    };

    SECTION("2d_cubic")           { test_quad_meshes<Symmetry::Cubic<>>(meshes, pattern_2d);           }
    // SECTION("2d_orthotropic")     { test_quad_meshes<Symmetry::Orthotropic<>>(meshes, pattern_2d);     }
    // SECTION("2d_diagonal")        { test_quad_meshes<Symmetry::Diagonal<>>(meshes, pattern_2d);        }
    // SECTION("2d_doubly_periodic") { test_quad_meshes<Symmetry::DoublyPeriodic<>>(meshes, pattern_2d);  }
    // SECTION("2d_square")          { test_quad_meshes<Symmetry::Square<>>(meshes, pattern_2d);          }
}

#endif
