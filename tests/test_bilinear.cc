////////////////////////////////////////////////////////////////////////////////
#include <isosurface_inflator/BilinearMap.hh>
#include <catch2/catch.hpp>
#include <random>
////////////////////////////////////////////////////////////////////////////////

void test_identity(const BilinearMap &f, Eigen::RowVector2d pts) {
    auto res = f.apply(pts[0], pts[1]);
    REQUIRE(std::abs(res[2]) < 1e-16);
    Eigen::RowVector2d q(res[0], res[1]);
    REQUIRE((q - pts).squaredNorm() < 1e-10);

    auto jac = f.jacobian(pts[0], pts[1]);
    REQUIRE((jac - Eigen::Matrix3d::Identity()).norm() < 1e-10);
}

////////////////////////////////////////////////////////////////////////////////

TEST_CASE("bilinear_map", "[bilinear_map]") {

    BilinearMap f;

    // Verify that the default bilinear map is indeed identity and that its
    // Jacobian is constant and equal to the identity.

    SECTION("fixed") {
        Eigen::Matrix<double, 5, 2> pts;
        pts <<
            -1, -1,
             1, -1,
             1,  1,
            -1,  1,
             0,  0;
        for (int i = 0; i < pts.rows(); ++i) {
            test_identity(f, pts.row(i));
        }
    }

    SECTION("random") {
        std::default_random_engine gen;
        std::uniform_real_distribution<double> dist(-10, 10);
        for (int i = 0; i < 100; ++i) {
            Eigen::RowVector2d p;
            p << dist(gen), dist(gen);
            test_identity(f, p);
        }
    }

}
