#include "RBFOrthoInflator.hh"

using namespace std;

RBFOrthoInflator::RBFOrthoInflator(Real epsilon, size_t dim) : m_dim(dim), m_rbfInflator(epsilon, dim + dim - 1) {
}

RBFOrthoInflator::RBFOrthoInflator(std::string png_path, Real epsilon, size_t dim) : m_dim(dim), m_rbfInflator(png_path, epsilon, dim + dim - 1) {
}

// VoxelsInflator methods
void
RBFOrthoInflator::m_inflate(const std::vector<Real> &reducedParams) {
    size_t totalDim = m_dim + m_dim - 1;

    std::vector<Real> params(totalDim * totalDim);
    for (size_t rp = 0; rp < reducedParams.size(); rp++) {
        vector<size_t> gluedParamIndices = reducedParamToAll(rp);

        for (unsigned j = 0; j < gluedParamIndices.size(); j++) {
            params[gluedParamIndices[j]] = reducedParams[rp];
        }
    }

    m_rbfInflator.inflate(params);

    m_vertices = m_rbfInflator.vertices();
    m_elements = m_rbfInflator.elements();
}

std::vector<VectorField<Real, 2>>
RBFOrthoInflator::volumeShapeVelocities() const {
    //std::cout << "[RBFOrthoInflator] Starting volumeShapeVelocities()" << std::endl;
    std::vector<VectorField<Real, N>> result;
    std::vector<VectorField<Real, N>> completeResult = m_rbfInflator.volumeShapeVelocities();
    size_t totalDim = m_dim + m_dim - 1;

    //std::cout << "[RBFOrthoInflator] After computing volumeShapeVelocities() from RBF" << std::endl;
    // For each parameter, we should add contribution for all related parameters in the volumeShapeVelocities
    for (size_t rp = 0; rp < m_dim * m_dim; rp++) {
        vector<size_t> gluedParamIndices = reducedParamToAll(rp);

        VectorField<Real, N> velocity = completeResult[gluedParamIndices[0]];
        for (unsigned j=1; j < gluedParamIndices.size(); j++) {
            size_t p = gluedParamIndices[j];
            velocity += completeResult[p];
        }

        result.push_back(velocity);
    }
    //std::cout << "[RBFOrthoInflator] Ending volumeShapeVelocities()" << std::endl;

    return result;
}

size_t
RBFOrthoInflator::numParameters() const {
    return m_dim * m_dim;
}

void
RBFOrthoInflator::loadMeshingOptions(const std::string &moptsPath) {
    m_rbfInflator.loadMeshingOptions(moptsPath);
}

void
RBFOrthoInflator::setMaxElementVolume(Real maxElementVol) {
    m_rbfInflator.setMaxElementVolume(maxElementVol);
}

Real
RBFOrthoInflator::getMaxElementVolume() const {
    return m_rbfInflator.getMaxElementVolume();
}

// Implementation of auxiliary functions
std::vector<std::vector<Real>>
RBFOrthoInflator::vecToMat(const std::vector<Real> &vec, size_t nRows, size_t nCols) const {
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
RBFOrthoInflator::matToVec(const std::vector<std::vector<Real>> &mat) const {
    std::vector<Real> result;

    for (size_t i = 0; i < mat.size(); i++) {
        result.insert(result.end(), mat[i].begin(), mat[i].end() );
    }

    return result;
}

std::vector<size_t>
RBFOrthoInflator::reducedParamToAll(size_t param) const {
    std::vector<size_t> result;
    size_t totalDim = m_dim + m_dim - 1;

    size_t ri = param / m_dim;
    size_t rj = param -  m_dim * ri;

    std::vector<size_t> is;
    std::vector<size_t> js;

    is.push_back(ri);
    if (ri < (m_dim - 1)) {
        is.push_back(totalDim - 1 - ri);
    }

    js.push_back(rj);
    if (rj < (m_dim - 1)) {
        js.push_back(totalDim - 1 - rj);
    }

    for (auto i : is)
        for (auto j : js)
            result.push_back(i * totalDim + j);

    return result;
}

size_t
RBFOrthoInflator::originalToReducedParam(size_t originalParam) const {
    size_t result;
    size_t totalDim = m_dim + m_dim - 1;

    size_t i = originalParam / totalDim;
    size_t j = originalParam - totalDim * i;

    size_t ri, rj;

    if (i > totalDim/2) {
        ri = totalDim - 1 - i;
    }
    else if (j > totalDim/2) {
        rj = totalDim - 1 - j;
    }

    result = ri * m_dim + rj;

    return result;
}

void
RBFOrthoInflator::savePng(const std::vector<Real> &reducedParams, std::string png_path) const {
    size_t totalDim = m_dim + m_dim - 1;

    std::vector<Real> params(totalDim * totalDim);
    for (size_t rp = 0; rp < reducedParams.size(); rp++) {
        vector<size_t> gluedParamIndices = reducedParamToAll(rp);

        for (unsigned j = 0; j < gluedParamIndices.size(); j++) {
            params[gluedParamIndices[j]] = reducedParams[rp];
        }
    }

    m_rbfInflator.savePng(params, png_path);
}
