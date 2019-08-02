#include "RBFInflator.hh"

using namespace std;

RBFInflator::RBFInflator(Real epsilon, size_t dim1, size_t dim2) {
    m_dim1 = dim1;
    m_dim2 = dim2;
    m_epsilon = epsilon;

    // dim1 relates to x, so number of columns
    // dim2 relates to y, so number of rows
    m_coeffMatrix.resize(dim2);
    for (size_t i = 0; i < dim2; i++) {
        m_coeffMatrix[i].resize(dim1);
    }
}

RBFInflator::RBFInflator(std::string png_path, Real epsilon, size_t dim1, size_t dim2) {
    m_dim1 = dim1;
    m_dim2 = dim2;
    m_epsilon = epsilon;

    RBF<Real> levelSet = RBF<Real>(png_path, epsilon, dim1, dim2);

    m_coeffMatrix = levelSet.coefficients();
}

// VoxelsInflator methods
void
RBFInflator::m_inflate(const std::vector<Real> &params) {
    // Makes sure we use the last meshing options set
    m_mesher.meshingOptions = m_meshingOptions;

    // Makes sure we use periodic mesher
    m_mesher.meshInterfaceConsistently = true;

    std::vector<std::vector<Real>> coeffMatrix = vecToMat(params, m_dim1, m_dim2);

    RBF<Real> levelSet = RBF<Real>(coeffMatrix, m_epsilon);

    m_mesher.meshSlice(levelSet, m_vertices, m_elements);
}

std::vector<VectorField<Real, 2>>
RBFInflator::volumeShapeVelocities() const {
    std::vector<VectorField<Real, 2>> result;

    throw std::runtime_error("TODO: implement volumeShapeVelocities() for RBF Inflator");

    return result;
}

size_t
RBFInflator::numParameters() const {
    return m_coeffMatrix.size() * m_coeffMatrix[0].size();
}

void
RBFInflator::loadMeshingOptions(const std::string &moptsPath) {
    m_meshingOptions.load(moptsPath);
}

void
RBFInflator::setMaxElementVolume(Real maxElementVol) {
    m_meshingOptions.maxArea = maxElementVol;
}

Real
RBFInflator::getMaxElementVolume() const {
    return m_meshingOptions.maxArea;
}


// Implementation of auxiliary functions
std::vector<std::vector<Real>>
RBFInflator::vecToMat(const std::vector<Real> &vec, size_t nRows, size_t nCols) const {
    std::vector<std::vector<Real>> result(nRows);

    for (size_t i = 0 ; i < nRows; i++) {
        vector<Real>::const_iterator first = vec.begin() + i * nCols;
        vector<Real>::const_iterator last = vec.begin() + (i+1) * nCols;
        vector<Real> newVec(first, last);
        result[i] = newVec;
    }

    return result;
}

std::vector<Real>
RBFInflator::matToVec(const std::vector<std::vector<Real>> &mat) const {
    std::vector<Real> result;

    for (size_t i = 0; i < mat.size(); i++) {
        result.insert(result.end(), mat[i].begin(), mat[i].end() );
    }

    return result;
}