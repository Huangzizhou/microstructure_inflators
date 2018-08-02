////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/IsosurfaceInflator.hh>
#include <isosurface_inflator/MeshingOptions.hh>
#include <isosurface_inflator/IsosurfaceInflatorConfig.hh>
#include <MeshFEM/StringUtils.hh>
#include <catch2/catch.hpp>
////////////////////////////////////////////////////////////////////////////////

void test_inflation(const std::string &mesher, const std::string &pattern) {
    size_t inflation_graph_radius = 2;

    IsosurfaceInflator inflator(mesher, true, pattern, inflation_graph_radius);
    std::vector<Real> params(inflator.defaultParameters());

    const nlohmann::json meshing_options = R"({
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
    inflator.meshingOptions().load(meshing_options);

    // inflator.disablePostprocess();
    // inflator.enableCheapPostprocess();
    // inflator.setReflectiveInflator(!args.nonReflectiveInflator);
    // inflator.setGenerateFullPeriodCell(!args.ortho_cell);

    std::cout << "Dilating " << mesher << ":\n\t" << pattern << std::endl;
    inflator.inflate(params);

    // if (MeshFEM::startswith(mesher, "2d_")) {
    //     MeshIO::save(mesher + ".obj", inflator.vertices(), inflator.elements());
    // } else {
    //     MeshIO::save(mesher + ".mesh", inflator.vertices(), inflator.elements());
    // }

    // Verify that normals have a zero z component in 2D
    if (MeshFEM::startswith(mesher, "2d_")) {
        const auto &n = inflator.vertexNormals();
        double maxZMag = 0;
        size_t maxZMagVtx = 0;
        for (size_t vi = 0; vi < n.size(); ++vi) {
            Real zmag = std::abs(n[vi][2]);
            if (zmag > maxZMag) {
                maxZMag = zmag;
                maxZMagVtx = vi;
            }
        }

        if (maxZMag > 0.1) {
            std::cerr << "Large normal z component: " << maxZMag << std::endl;
            auto n = inflator.trackSignedDistanceGradient(inflator.vertices().at(maxZMagVtx).point);
            std::cout << "normal: " << n.transpose() << std::endl;
        }
        REQUIRE(maxZMag <= 0.1);
    }
}

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("inflate_default_pattern", "[isosurface_inflation]") {
    auto pattern_2d = DATA_DIR "patterns/2D/topologies/0001.obj";
    auto pattern_3d = DATA_DIR "patterns/3D/reference_wires/pattern0000.wire";

    SECTION("2d_cubic")           { test_inflation("2d_cubic", pattern_2d);           }
    SECTION("2d_orthotropic")     { test_inflation("2d_orthotropic", pattern_2d);     }
    SECTION("2d_diagonal")        { test_inflation("2d_diagonal", pattern_2d);        }
    // SECTION("2d_non_periodic")    { test_inflation("2d_non_periodic", pattern_2d);    } // Not working
    SECTION("2d_doubly_periodic") { test_inflation("2d_doubly_periodic", pattern_2d); }
    SECTION("2d_square")          { test_inflation("2d_square", pattern_2d);          }
    SECTION("3d_cubic")           { test_inflation("cubic", pattern_3d);              }
    SECTION("3d_orthotropic")     { test_inflation("orthotropic", pattern_3d);        }
    // SECTION("3d_non_periodic")    { test_inflation("non_periodic", pattern_3d);       } // Not working
    // SECTION("3d_triply_periodic") { test_inflation("triply_periodic", pattern_3d);    } // Needs CGAL triply periodic meshing (4.13)
    SECTION("3d_square")          { test_inflation("square", pattern_3d);             }
}
