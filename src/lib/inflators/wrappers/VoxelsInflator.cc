#include "VoxelsInflator.hh"

using namespace std;

VoxelsInflator::VoxelsInflator(const std::vector<std::vector<Real>> &densityMatrix) : m_densityMatrix(densityMatrix) {
    m_nRows = densityMatrix.size();
    m_nCols = densityMatrix[0].size();
}


// VoxelsInflator methods
void
VoxelsInflator::m_inflate(const std::vector<Real> &params) {
    // Makes sure we use the last meshing options set
    m_mesher.meshingOptions = m_meshingOptions;

    // Makes sure we use periodic mesher
    m_mesher.meshInterfaceConsistently = true;

    std::vector<std::vector<Real>> densityMatrix = vecToMat(params, m_nRows, m_nCols);

    VoxelSDF<Real> sdf = VoxelSDF<Real>(densityMatrix);

    m_mesher.meshSlice(sdf, m_vertices, m_elements);
}

std::vector<VectorField<Real, 2>>
VoxelsInflator::volumeShapeVelocities() const {
    std::vector<VectorField<Real, 2>> result;

    // Does it make sense though?
    throw std::runtime_error("TODO: implement volumeShapeVelocities() for Voxels Inflator");

    return result;
}

size_t
VoxelsInflator::numParameters() const {
    return m_densityMatrix.size() * m_densityMatrix[0].size();
}

void
VoxelsInflator::loadMeshingOptions(const std::string &moptsPath) {
    m_meshingOptions.load(moptsPath);
}

void
VoxelsInflator::setMaxElementVolume(Real maxElementVol) {
    m_meshingOptions.maxArea = maxElementVol;
}

Real
VoxelsInflator::getMaxElementVolume() const {
    return m_meshingOptions.maxArea;
}


// Implementation of auxiliary functions
std::vector<std::vector<Real>>
VoxelsInflator::vecToMat(const std::vector<Real> &vec, size_t nRows, size_t nCols) {
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
VoxelsInflator::matToVec(const std::vector<std::vector<Real>> &mat) {
    std::vector<Real> result;

    for (size_t i = 0; i < mat.size(); i++) {
        result.insert(result.end(), mat[i].begin(), mat[i].end() );
    }

    return result;
}