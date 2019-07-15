////////////////////////////////////////////////////////////////////////////////
// VoxelSDF.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computes signed distance represented by voxel densities
*/
////////////////////////////////////////////////////////////////////////////////
#ifndef VOXELSIGNEDDISTANCE_HH
#define VOXELSIGNEDDISTANCE_HH

#include "../isosurface_inflator/SignedDistanceRegion.hh"
#include "../isosurface_inflator/SignedDistance.hh"

#include <MeshFEM/Future.hh>
#include <unordered_map>
#include <unordered_set>
#include <math.h>


namespace SD = SignedDistance;

template<typename _Real>
class VoxelSDF : public SignedDistanceRegion<2> {

public:

    using Real = _Real;
    VoxelSDF(const std::vector<std::vector<Real>> &densityMatrix);

    // Always support double type for compatibility with SignedDistanceRegion
    virtual double signedDistance(const Point2D &p) const override {
        return stripAutoDiff(m_signedDistanceImpl(autodiffCast<Point2<Real>>(p)));
    }

    // Also support automatic differentiation types
    template<typename Real2, bool DebugDerivatives = false>
    Real2 signedDistance(const Point2<Real2> &p) const { return m_signedDistanceImpl(p); }

    size_t numParams() const { return m_densityMatrix.size() * m_densityMatrix[0].size(); }

    // Representative cell bounding box (region to be meshed)
    virtual const BBox<Point2D> &boundingBox() const override { return m_bbox; }

    virtual ~VoxelSDF() override = default;


private:

    // Additional Real type to support automatic differentiation wrt. p only
    template<typename Real2>
    Real2 m_signedDistanceImpl(const Point2<Real2> p) const;

    // Bounding box for the meshing cell.
    BBox<Point2D> m_bbox;

    // density matrix
    std::vector<std::vector<Real>> m_densityMatrix;
};

#endif /* end of include guard: VOXELSIGNEDDISTANCE_HH */
