//
// Created by Davi Colli Tozoni on 3/26/18.
//

#include "ConstrainedIsoinflator.hh"

using namespace std;

template<size_t N>
ConstrainedIsoinflator<N>::ConstrainedIsoinflator(const std::string &wireMeshPath, const std::string &symmetryType, bool vertex_thickness, const std::vector<bool> &paramsMask, const std::vector<double> &originalParams, size_t inflationGraphRadius) {
    m_originalParams = originalParams;

    // Creates map from constrained parameters to original ones
    size_t originalIdx = 0;
    m_numFilteredParams = 0;
    for (; originalIdx < paramsMask.size(); originalIdx++) {
        if (!paramsMask[originalIdx]) {
            m_filteredParamsToOriginal.push_back(originalIdx);
            m_numFilteredParams++;
        }
    }

    // Create and save inflator
    m_infl = Future::make_unique<IsoinflatorWrapper<N>>(wireMeshPath, symmetryType, vertex_thickness, inflationGraphRadius);
}

template<size_t N>
void ConstrainedIsoinflator<N>::m_inflate(const std::vector<Real> &params) {
    m_infl->inflate(constrainedToFullParameters(params));
}

template<size_t N>
bool ConstrainedIsoinflator<N>::isPrintable(const std::vector<Real> &params) {
    return m_infl->isPrintable(constrainedToFullParameters(params));
}

template<size_t N>
std::vector<VectorField<Real, N>> ConstrainedIsoinflator<N>::volumeShapeVelocities() const {
    std::vector<VectorField<Real, N>> result;
    std::vector<VectorField<Real, N>> originalVelocities = m_infl->volumeShapeVelocities();

    for (unsigned i=0; i < numParameters(); i++) {
        result.push_back(originalVelocities[m_filteredParamsToOriginal[i]]);
    }

    return result;
}

template<size_t N>
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
ConstrainedIsoinflator<N>::selfSupportingConstraints(const std::vector<double> &params) const {
    return m_infl->selfSupportingConstraints(constrainedToFullParameters(params));
}

template<size_t N>
size_t ConstrainedIsoinflator<N>::numParameters() const {
    return m_numFilteredParams;
}

template<size_t N>
ParameterType ConstrainedIsoinflator<N>::parameterType(size_t p) const {
    return m_infl->parameterType(m_filteredParamsToOriginal[p]);
}

template<size_t N>
template<typename T>
std::vector<T> ConstrainedIsoinflator<N>::constrainedToFullParameters(const std::vector<T> &constrainedParameters) const {
    vector<T> result(m_originalParams.size());

    // copy original vector
    for (unsigned i = 0; i < m_originalParams.size(); i++) {
        result[i] = m_originalParams[i];
    }

    // change only the constrained parameters using our map
    for (unsigned i=0; i < numParameters(); i++) {
        result[m_filteredParamsToOriginal[i]] = constrainedParameters[i];
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Explicit Instantiations: 2D and 3D inflators.
////////////////////////////////////////////////////////////////////////////////
template class ConstrainedIsoinflator<2>;
template class ConstrainedIsoinflator<3>;