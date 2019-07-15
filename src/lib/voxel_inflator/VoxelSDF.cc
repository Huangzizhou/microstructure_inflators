#include "VoxelSDF.hh"
#include <unsupported/Eigen/AutoDiff>

using PVec = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using ADScalar = Eigen::AutoDiffScalar<PVec>;

template<typename Real>
VoxelSDF<Real>::VoxelSDF(const std::vector<std::vector<Real>> &densityMatrix) : m_densityMatrix(densityMatrix) {
    Point2D min = Point2D(-1.0, -1.0);
    Point2D max = Point2D( 1.0,  1.0);

    m_bbox = BBox<Point2D>(min, max);
}


// Remember: matrix is defined row by row. So, row index corresponds to y and column index to x
template<typename Real>
template<typename Real2>
Real2 VoxelSDF<Real>::m_signedDistanceImpl(const Point2<Real2> p) const {
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

////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations: double (need some more implementation for auto diff)
////////////////////////////////////////////////////////////////////////////////
template class VoxelSDF<double>;
//template class VoxelSDF<ADScalar>;
