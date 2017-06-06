#ifndef IGLSURFACEMESHERMC_HH
#define IGLSURFACEMESHERMC_HH

#include "MesherBase.hh"
#include <MeshIO.hh>

#include <igl/copyleft/marching_cubes.h>

#include <Parallelism.hh>

class IGLSurfaceMesherMC : public MesherBase {
public:
    using Real = SignedDistanceRegion<3>::Real;
    using MesherBase::MesherBase;
    using MesherBase::meshingOptions;

    virtual void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements) override {
        mesh(sdf, vertices, elements, 0.0);
    }

    void mesh(const SignedDistanceRegion<3> &sdf,
            std::vector<MeshIO::IOVertex> &vertices,
            std::vector<MeshIO::IOElement> &elements,
            const double isolevel);
};

mesh(const SignedDistanceRegion<3> &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &elements, const double isolevel)
{
    size_t gs = meshingOptions.marchingCubesGridSize;
    const auto bbox = sdf.boundingBox();

    const size_t nsamples = gs * gs * gs;

    // Evaluation point locations;
    // flattened to be accessed as:
    // xi + gs * (yi + gs * zi)
    Eigen::MatrixXd sampleLocations(nsamples, 3);
    {
        size_t i = 0;
        for (size_t zi = 0; zi < gs; ++zi) {
            for (size_t yi = 0; yi < gs; ++yi) {
                for (size_t xi = 0; xi < gs; ++xi) {
                    sampleLocations.row(i) = bbox.interpolatePoint(
                            Point3D(xi / Real(gs - 1.0),
                                    yi / Real(gs - 1.0),
                                    zi / Real(gs - 1.0)));
                    ++i;
                }
            }
        }
    }

    // Evaluate signed distances at each grid point
    Eigen::VectorXd signedDistances(nsamples);
#if USE_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nsamples),
            [&](const tbb::blocked_range<size_t> &r) {
                for (size_t i = r.begin(); i < r.end(); ++i)
                    signedDistances(i) = sdf.signedDistance(sampleLocations.row(i)) - isolevel;
            }
        );
#else
    for (size_t i = 0; i < nsamples; ++i)
        signedDistances(i) = sdf.signedDistance(sampleLocations.row(i)) - isolevel;
#endif
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::copyleft::marching_cubes(signedDistances, sampleLocations,
            gs, gs, gs, V, F);

    vertices.clear(), elements.clear();
    vertices.reserve(V.rows());
    elements.reserve(F.rows());
    for (size_t i = 0; i < V.rows(); ++i)
        vertices.emplace_back(Point3D(V.row(i)));
    for (size_t i = 0; i < F.rows(); ++i)
        elements.emplace_back(F(i, 0), F(i, 1), F(i, 2));
}

#endif /* end of include guard: IGLSURFACEMESHERMC_HH */
