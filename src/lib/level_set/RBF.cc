#define STB_IMAGE_IMPLEMENTATION

#include "RBF.hh"
#include <unsupported/Eigen/AutoDiff>

using PVec = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
using ADScalar = Eigen::AutoDiffScalar<PVec>;

template<typename Real>
RBF<Real>::RBF(const std::vector<std::vector<Real>> &coeffMatrix, Real epsilon) : m_coeffMatrix(coeffMatrix) {
    Point2D min = Point2D(-1.0, -1.0);
    Point2D max = Point2D( 1.0,  1.0);

    m_bbox = BBox<Point2D>(min, max);

    m_d1 = m_coeffMatrix.size();
    m_d2 = m_coeffMatrix[0].size();

    m_min1 = -1.0;
    m_max1 =  1.0;
    m_min2 = -1.0;
    m_max2 =  1.0;

    m_epsilon = epsilon;

    m_dt1 = (m_max1 - m_min1) / (m_d1 - 1);
    m_dt2 = (m_max2 - m_min2) / (m_d2 - 1);
}

////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations: double (need some more implementation for auto diff)
////////////////////////////////////////////////////////////////////////////////
template class RBF<double>;
template class RBF<ADScalar>;
