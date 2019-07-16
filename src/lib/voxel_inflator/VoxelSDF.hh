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

    // Remember: matrix is defined row by row. So, row index corresponds to y and column index to x
    // Additional Real type to support automatic differentiation wrt. p only
    template<typename Real2>
    Real2 m_signedDistanceImpl(const Point2<Real2> p) const {
        size_t cells_per_axis = m_densityMatrix.size();

        // Base cell is from -1 to 1.
        Real2 x = p[0]; // on x axis
        Real2 y = p[1]; // on y axis

        size_t i = floor(cells_per_axis * (y + 1.0) / 2.0);
        size_t j = floor(cells_per_axis * (x + 1.0) / 2.0);

        if (y == 1.0)
            i = m_densityMatrix.size() - 1;

        if (x == 1.0)
            j = m_densityMatrix[0].size() - 1;

        m_densityMatrix[i][j];
        Real2 density = m_densityMatrix[i][j];

        // Mapping density from 0..1 to 1..-1
        Real2 distance = -1.0 * ((density * 2.0) - 1.0);

        return distance;
    }

    // Bounding box for the meshing cell.
    BBox<Point2D> m_bbox;

    // density matrix
    std::vector<std::vector<Real>> m_densityMatrix;
};

#endif /* end of include guard: VOXELSIGNEDDISTANCE_HH */
