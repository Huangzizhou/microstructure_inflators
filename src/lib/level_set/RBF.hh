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

    RBF(std::string png_path, Real epsilon, size_t dim1, size_t dim2) {
        Point2D min = Point2D(-1.0, -1.0);
        Point2D max = Point2D( 1.0,  1.0);

        m_bbox = BBox<Point2D>(min, max);

        m_d1 = dim1;
        m_d2 = dim2;

        m_min1 = -1.0;
        m_max1 =  1.0;
        m_min2 = -1.0;
        m_max2 =  1.0;

        m_epsilon = epsilon;

        m_dt1 = (m_max1 - m_min1) / (m_d1 - 1);
        m_dt2 = (m_max2 - m_min2) / (m_d2 - 1);

        int width, height;
        int channels;
        unsigned char * png_matrix = stbi_load(png_path.c_str(), &width, &height, &channels, STBI_grey);

        // Construct density matrix
        std::vector<Real> data_values(height * width);

        if (channels == 1) {
            for (size_t y = 0; y < height; y++) {
                unsigned char *row = &png_matrix[y * width];

                for (size_t x = 0; x < width; x++) {
                    unsigned char byte = row[x];
                    data_values[y*width + x] = byte / 255.0 * 2.0 - 1.0;
                    std::cout << (data_values[y*width + x] + 1) / 2;
                }
                std::cout << std::endl;
            }
        }

        size_t num_rows = height;
        size_t num_cols = width;

        // Construct grid with positions between -1 to 1 on both directions
        // Given pixel size, first data point should be at -1 + pixel_size/2. Last should be at 1 - pixel_size/2. All the other N - 2
        // in between
        Real pixel_size1 = (1.0 - (-1.0)) / num_cols;
        Real pixel_size2 = (1.0 - (-1.0)) / num_rows;
        std::vector<Real> x1;
        std::vector<Real> x2;
        for (size_t i = 0; i < num_rows; i++) {
            Real x2_pos = (1.0 - pixel_size2 / 2.0) - i * pixel_size2;
            for (size_t j = 0; j < num_cols; j++) {
                Real x1_pos = (-1.0 + pixel_size1 / 2.0) + j * pixel_size1;

                x1.push_back(x1_pos);
                x2.push_back(x2_pos);
            }
        }

        std::cout << "Running fitting..." << std::endl;
        std::vector<Real> coefficients = fit<Real>(x1, x2, data_values, dim1, dim2, epsilon);
        std::cout << "Ending fitting..." << std::endl;

        m_coeffMatrix.resize(dim2);
        for (size_t i = 0 ; i < dim2; i++) {
            typename std::vector<Real>::const_iterator first = coefficients.begin() + i * dim1;
            typename std::vector<Real>::const_iterator last = coefficients.begin() + (i+1) * dim1;
            std::vector<Real> newVec(first, last);
            m_coeffMatrix[i] = newVec;
        }
    }

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
