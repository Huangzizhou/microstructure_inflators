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

template<typename Real>
RBF<Real>::RBF(std::string png_path, Real epsilon, size_t dim1, size_t dim2) {
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
        for (size_t y = 0; y < (size_t) height; y++) {
            unsigned char *row = &png_matrix[y * width];

            for (size_t x = 0; x < (size_t) width; x++) {
                unsigned char byte = row[x];
                data_values[y*width + x] = byte / 255.0 * 2.0 - 1.0;
                //std::cout << (data_values[y*width + x] + 1) / 2;
            }
            //std::cout << std::endl;
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

    //std::cout << "Running fitting..." << std::endl;
    std::vector<Real> coefficients = fit<Real>(x1, x2, data_values, dim1, dim2, epsilon);
    //std::cout << "Ending fitting..." << std::endl;

    m_coeffMatrix.resize(dim2);
    for (size_t i = 0 ; i < dim2; i++) {
        typename std::vector<Real>::const_iterator first = coefficients.begin() + i * dim1;
        typename std::vector<Real>::const_iterator last = coefficients.begin() + (i+1) * dim1;
        std::vector<Real> newVec(first, last);
        m_coeffMatrix[i] = newVec;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations: double (need some more implementation for auto diff)
////////////////////////////////////////////////////////////////////////////////
template class RBF<double>;
template class RBF<ADScalar>;
