#include "VCGSurfaceMesher.hh"

#include "WireMesh.hh"
#include "PatternSignedDistance.hh"

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/create/marching_cubes.h>
#include <vcg/complex/algorithms/create/extended_marching_cubes.h>
#include <vcg/complex/algorithms/create/mc_trivial_walker.h>

class MyFace;
class MyVertex;

struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>::AsVertexType, vcg::Use<MyFace>::AsFaceType>{};
class  MyVertex    : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags>{};
class  MyFace      : public vcg::Face< MyUsedTypes, vcg::face::VertexRef, vcg::face::BitFlags> {};
class  MyMesh      : public vcg::tri::TriMesh< std::vector< MyVertex>, std::vector< MyFace > > {};

template<class SignedDistanceFunction>
void VCGSurfaceMesher<SignedDistanceFunction>::
mesh(const SignedDistanceFunction &sdf,
     std::vector<MeshIO::IOVertex> &vertices,
     std::vector<MeshIO::IOElement> &elements)
{
    typedef vcg::SimpleVolume<vcg::SimpleVoxel<Real>> MyVolume;
    typedef vcg::tri::TrivialWalker<MyMesh, MyVolume>    MyWalker;
    typedef vcg::tri::MarchingCubes<MyMesh, MyWalker>    MyMarchingCubes;

    size_t gs = meshingOptions.marchingCubesGridSize;
    auto bbox = SignedDistanceFunction::boundingBox();
    auto bbmin = bbox.minCorner, bbmax = bbox.maxCorner;

    // Ugly hack until VCG fixes its off-by-one (two!!??) error
    bbmax *= (gs + 2.0) / gs;
    vcg::Box3d bb(vcg::Point3d(bbmin[0], bbmin[1], bbmin[2]),
                  vcg::Point3d(bbmax[0], bbmax[1], bbmax[2]));
    
    MyVolume volume;
    volume.Init(vcg::Point3i(gs, gs, gs), bb);
    for (size_t i = 0; i < gs; ++i) {
        for (size_t j = 0; j < gs; ++j) {
            for (size_t k = 0; k < gs; ++k) {
                // Ugly hack until VCG fixes its off-by-one (two!!??) error
                auto p = bbox.interpolatePoint(Point3D(i / Real(gs - 2),
                                                       j / Real(gs - 2),
                                                       k / Real(gs - 2)));

                volume.Val(i, j, k) = sdf.signedDistance(p);
            }
        }
    }

    // MARCHING CUBES
    MyMesh mc_mesh;
    MyWalker walker;
    MyMarchingCubes mc(mc_mesh, walker);
    walker.template BuildMesh<MyMarchingCubes>(mc_mesh, volume, mc, 0);
    // Note: ignores the DELETED vertex flag (accessed through IsD())
    for (const auto &v : mc_mesh.vert)
        vertices.emplace_back(v.P()[0], v.P()[1], v.P()[2]);
    for (      auto &f : mc_mesh.face) { // f.V() isn't const :(
        if (f.VN() != 3) throw std::runtime_error("Expected triangles from marching cubes");
        elements.emplace_back(vcg::tri::Index(mc_mesh, f.V(0)),
                              vcg::tri::Index(mc_mesh, f.V(1)),
                              vcg::tri::Index(mc_mesh, f.V(2)));
    }
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiations
////////////////////////////////////////////////////////////////////////////////
template class VCGSurfaceMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Cubic<>>>>;
// Enable for slower builds...
// template class VCGSurfaceMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>>>>;
// template class VCGSurfaceMesher<PatternSignedDistance<double, WireMesh<ThicknessType::Vertex, Symmetry::TriplyPeriodic<>>>>;
