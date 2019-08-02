////////////////////////////////////////////////////////////////////////////////
// RBF.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computes signed distance represented by voxel densities
*/
////////////////////////////////////////////////////////////////////////////////
#ifndef RBF_HH
#define RBF_HH

#include <isosurface_inflator/SignedDistanceRegion.hh>
#include <isosurface_inflator/SignedDistance.hh>

#include <MeshFEM/Future.hh>
#include <unordered_map>
#include <unordered_set>
#include <math.h>

#include "../../3rdparty/libigl/external/stb/stb_image.h"

namespace SD = SignedDistance;

template<typename _Real>
class RBF : public SignedDistanceRegion<2> {

public:

    using Real = _Real;
    RBF(const std::vector<std::vector<Real>> &coeffMatrix, Real epsilon);

    RBF(std::string png_path, Real epsilon, size_t dim1, size_t dim2);

    // Always support double type for compatibility with SignedDistanceRegion
    virtual double signedDistance(const Point2D &p) const override {
        return stripAutoDiff(m_signedDistanceImpl(autodiffCast<Point2<Real>>(p)));
    }

    // Also support automatic differentiation types
    template<typename Real2, bool DebugDerivatives = false>
    Real2 signedDistance(const Point2<Real2> &p) const { return m_signedDistanceImpl(p); }

    size_t numParams() const { return m_coeffMatrix.size() * m_coeffMatrix[0].size(); }

    // Representative cell bounding box (region to be meshed)
    virtual const BBox<Point2D> &boundingBox() const override { return m_bbox; }

    virtual ~RBF() override = default;

    std::vector<std::vector<Real>> coefficients() {
        return m_coeffMatrix;
    }


private:

    // Additional Real type to support automatic differentiation wrt. p only
    template<typename Real2>
    Real2 m_signedDistanceImpl(const Point2<Real2> p) const {
        Real2 result = 0.0;

        for (size_t i=0; i<m_d1; i++) {
            for (size_t j=0; j<m_d2; j++) {
                // Find rbf point correspond to i and j
                Real xi = m_min1 + i * m_dt1;
                Real xj = m_min2 + j * m_dt1;

                Real2 r = sqrt((p[0] - xi) * (p[0] - xi) + (p[1] - xj) * (p[1] - xj));

                result += m_coeffMatrix[i][j] * exp(- (m_epsilon * r) * (m_epsilon * r));
            }
        }

        return result;
    }

    template<typename Real2>
    std::vector<Real> fit(std::vector<Real> x1, std::vector<Real> x2, std::vector<Real> values, size_t d1 = 5, size_t d2 = 5, Real epsilon = 0.5) {
        // create rhs vector
        Eigen::Matrix<Real, Eigen::Dynamic, 1> B = Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>>(values.data(), values.size());

        // create matrix
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> A(x1.size(), d1 * d2);
        size_t data_size = x1.size();

        for (size_t ip = 0; ip < data_size; ip++) {
            for (size_t i = 0; i < d1; i++) {
                for (size_t j = 0; j < d2; j++) {
                    // Find rbf point correspond to i and j
                    Real xi = m_min1 + i * m_dt1;
                    Real xj = m_min2 + j * m_dt1;

                    Real2 r = sqrt((x1[ip] - xi) * (x1[ip] - xi) + (x2[ip] - xj) * (x2[ip] - xj));

                    A(ip, m_d2 * i + j) = exp(- (m_epsilon * r) * (m_epsilon * r));
                }
            }
        }

        // solving least squares
        //std::cout << "Solving" << std::endl;
        //Eigen::Matrix<Real, Eigen::Dynamic, 1> x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B);
        Eigen::Matrix<Real, Eigen::Dynamic, 1> x = A.colPivHouseholderQr().solve(B);
        //Eigen::Matrix<Real, Eigen::Dynamic, 1> x = (A.transpose() * A).ldlt().solve(A.transpose() * B);
        std::vector<Real> result(x.data(), x.data() + x.rows() * x.cols());

        return result;
    }

    // Bounding box for the meshing cell.
    BBox<Point2D> m_bbox;

    // density matrix
    std::vector<std::vector<Real>> m_coeffMatrix;
    size_t m_d1, m_d2;
    Real m_epsilon;
    Real m_dt1, m_dt2;
    Real m_min1, m_min2;
    Real m_max1, m_max2;
};

#endif /* end of include guard: VOXELSIGNEDDISTANCE_HH */
