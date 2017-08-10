#ifndef IGLSURFACEMESHERMC_HH
#define IGLSURFACEMESHERMC_HH

#include "MesherBase.hh"
#include <MeshIO.hh>

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

#endif /* end of include guard: IGLSURFACEMESHERMC_HH */
