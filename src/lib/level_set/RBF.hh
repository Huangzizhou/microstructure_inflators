////////////////////////////////////////////////////////////////////////////////
// RBF.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computes signed distance represented by voxel densities
*/
////////////////////////////////////////////////////////////////////////////////
#ifndef RBF_HH
#define RBF_HH

#define  DEBUG_OUT 1

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
    using Vector2D = Eigen::Matrix<Real, 2, 1>;
    RBF(const std::vector<std::vector<Real>> &coeffMatrix, Real epsilon);

    RBF(std::string png_path, Real epsilon, size_t dim);

    // Always support double type for compatibility with SignedDistanceRegion
    virtual double signedDistance(const Point2D &p) const override {
        return stripAutoDiff(m_signedDistanceImpl(autodiffCast<Point2<Real>>(p), m_coeffMatrix));
    }

    // Also support automatic differentiation types
    template<typename Real2, bool DebugDerivatives = false>
    Real2 signedDistance(const Point2<Real2> &p) const { return m_signedDistanceImpl(p, m_coeffMatrix); }

    size_t numParams() const { return m_coeffMatrix.size() * m_coeffMatrix[0].size(); }

    // Representative cell bounding box (region to be meshed)
    virtual const BBox<Point2D> &boundingBox() const override { return m_bbox; }

    virtual ~RBF() override = default;

    std::vector<std::vector<Real>> coefficients() {
        return m_coeffMatrix;
    }

    Real partialDerivative(size_t i, size_t j, Vector2D x) const {
        Real result = basis(i, j, x);

#if DEBUG_OUT
        Real currentSD = signedDistance(x);

        // perturb i j coefficient
        Real pert = 1e-6;
        std::vector<std::vector<Real>> perturbedCoeffMatrix = m_coeffMatrix;
        perturbedCoeffMatrix[i][j] += pert;
        Real perturbedSD = m_signedDistanceImpl(x, perturbedCoeffMatrix);
        Real finiteDiff = (perturbedSD - currentSD) / pert;

        if (abs(finiteDiff - result) > 1e-6) {
            std::cout << "[Warning!] x: " << x << std::endl;
            std::cout << "Partial derivative / finite diff: " << result << " / " << finiteDiff << std::endl;
            std::cout << "Error: " << abs(finiteDiff - result) << std::endl;
        }

        //std::cout << "Partial derivative / finite diff: " << result << " / " << finiteDiff << std::endl;
#endif

        return basis(i, j, x);
    }

    Vector2D gradient(Vector2D x) const {
        Vector2D result;
        result.setZero();

        for (size_t i = 0; i < m_dim; i++) {
            for (size_t j = 0; j < m_dim; j++) {
                Real xi = m_min + i * m_dt;
                Real xj = m_min + j * m_dt;

                Real r = sqrt((x(0) - xi) * (x(0) - xi) + (x(1) - xj) * (x(1) - xj));
                Vector2D diff;
                Real diff_x = x(0) - xi;
                Real diff_y = x(1) - xj;
                diff << diff_x, diff_y;

                Real w_ij = m_coeffMatrix[i][j];

                result += -2.0 * w_ij * m_epsilon * m_epsilon * exp(- (m_epsilon * r) * (m_epsilon * r)) * diff;
            }
        }

#if DEBUG_OUT
        Real currentSD = signedDistance(x);
        Real pert = 1e-6;

        // perturb x first coordinate
        Vector2D xPert0, xPert1;
        xPert0 << x(0) + pert, x(1);
        xPert1 << x(0), x(1) + pert;

        Real perturbedSD0 = signedDistance(xPert0);
        Real perturbedSD1 = signedDistance(xPert1);
        Vector2D finiteDiff;
        finiteDiff << (perturbedSD0 - currentSD) / pert, (perturbedSD1 - currentSD) / pert;

        if ((finiteDiff - result).norm() > 1e-5) {
            std::cout << "[Warning!] x: " << x << std::endl;
            std::cout << "[Warning!] Gradient / finite diff: " << result << " / " << finiteDiff << std::endl;
            std::cout << "Error: " << (finiteDiff - result).norm() << std::endl;
        }
#endif

        return result;
    }


private:

    Real basis(size_t i, size_t j, Vector2D x) const {
        Real xi = m_min + i * m_dt;
        Real xj = m_min + j * m_dt;

        Real r = sqrt((x(0) - xi) * (x(0) - xi) + (x(1) - xj) * (x(1) - xj));

        return exp(- (m_epsilon * r) * (m_epsilon * r));
    }

    // Additional Real type to support automatic differentiation wrt. p only
    template<typename Real2>
    Real2 m_signedDistanceImpl(const Point2<Real2> p, std::vector<std::vector<Real>> coeffMatrix) const {
        Real2 result = 0.0;

        for (size_t i=0; i<m_dim; i++) {
            for (size_t j=0; j<m_dim; j++) {
                // Find rbf point correspond to i and j
                Real xi = m_min + i * m_dt;
                Real xj = m_min + j * m_dt;

                Real2 r = sqrt((p[0] - xi) * (p[0] - xi) + (p[1] - xj) * (p[1] - xj));

                result += coeffMatrix[i][j] * exp(- (m_epsilon * r) * (m_epsilon * r));
            }
        }

        return result;
    }

    template<typename Real2>
    std::vector<Real> fit(std::vector<Real> x1, std::vector<Real> x2, std::vector<Real> values, size_t d = 5, Real epsilon = 0.5) {
        // create rhs vector
        Eigen::Matrix<Real, Eigen::Dynamic, 1> B = Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, 1>>(values.data(), values.size());

        // create matrix
        Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> A(x1.size(), d * d);
        size_t data_size = x1.size();

        for (size_t ip = 0; ip < data_size; ip++) {
            for (size_t i = 0; i < d; i++) {
                for (size_t j = 0; j < d; j++) {
                    // Find rbf point correspond to i and j
                    Real xi = m_min + i * m_dt;
                    Real xj = m_min + j * m_dt;

                    Real2 r = sqrt((x1[ip] - xi) * (x1[ip] - xi) + (x2[ip] - xj) * (x2[ip] - xj));

                    A(ip, d * i + j) = exp(- (m_epsilon * r) * (m_epsilon * r));
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
    size_t m_dim;
    Real m_epsilon;
    Real m_dt;
    Real m_min;
    Real m_max;
};

#endif /* end of include guard: VOXELSIGNEDDISTANCE_HH */
