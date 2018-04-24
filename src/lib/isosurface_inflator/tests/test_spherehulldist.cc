#include "SphereConvexHull.hh"
#include <MeshFEM/filters/gen_grid.hh>
#include <MeshFEM/GlobalBenchmark.hh>

using namespace SD::Primitives;

int main(int argc, const char *argv[]) {
    std::vector<Point3D> centers = {
        Point3D(0, 0, 0),
        Point3D(0.5, 0, 0),
        Point3D(0, 0.5, 0),
        Point3D(0, 0, 0.5) };

    std::vector<double> radii = { 0.1, 0.15, 0.25, 0.2 };

    SphereConvexHull<double> shull(centers, radii);


    {
        std::vector<MeshIO::IOVertex > gv;
        std::vector<MeshIO::IOElement> ge;
        size_t gs = 100;
        std::vector<size_t> grid_sizes({gs, gs, gs});
        gen_grid(grid_sizes, gv, ge);
        ScalarField<Real> sd(gv.size());
        BENCHMARK_START_TIMER("SD Sample");
        for (size_t i = 0; i < gv.size(); ++i) {
            for (size_t j = 0; j < 3; ++j) {
                gv[i].point[j] *= 2.0 / gs;
                gv[i].point[j] -= 1.0;
            }
            sd[i] = shull.signedDistance(gv[i].point);
        }
        BENCHMARK_STOP_TIMER("SD Sample");
        MSHFieldWriter writer("sd_debug.msh", gv, ge);
        writer.addField("sd", sd, DomainType::PER_NODE);
    }

    BENCHMARK_START_TIMER("CGAL SD Sample");
    shull.writeGroundTruth(256000, 100, "sd_cgal.msh");
    BENCHMARK_STOP_TIMER("CGAL SD Sample");

    shull.writeDebugInfo();

    BENCHMARK_REPORT();
    return 0;
}
