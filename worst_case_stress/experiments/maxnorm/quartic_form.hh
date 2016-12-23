////////////////////////////////////////////////////////////////////////////////
// quartic_form.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Evaluation of the quartic form (and its gradients/Hessians) represented
//      by a rank 4 tensor with major symmetries.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/17/2016 16:33:45
////////////////////////////////////////////////////////////////////////////////
#ifndef QUARTIC_FORM_HH
#define QUARTIC_FORM_HH

#include <ElasticityTensor.hh>
#include <SymmetricMatrix.hh>
#include <Types.hh>

template<typename Real, size_t N>
class QuarticForm {
public:
    QuarticForm(ElasticityTensor<Real, N> &T) : m_T(T) { }

    Real operator()(const VectorND<N> &x) const {
        Real result = 0;
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j)
                for (size_t k = 0; k < N; ++k)
                    for (size_t l = 0; l < N; ++l)
                        result += m_T(i, j, k, l) * x[i] * x[j] * x[k] * x[l];
        return result;
    }

    VectorND<N> gradient(const VectorND<N> &x) const {
        VectorND<N> g;
        g.setZero();
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j)
                for (size_t k = 0; k < N; ++k)
                    for (size_t l = 0; l < N; ++l)
                        g[i] += m_T(i, j, k, l) * x[j] * x[k] * x[l];
        return g;
    }

    SymmetricMatrixValue<Real, N> hessian(const VectorND<N> &x) const {
        SymmetricMatrixValue<Real, N> h;
        h.clear();
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < N; ++j)
                for (size_t k = 0; k < N; ++k)
                    for (size_t l = 0; l < N; ++l)
                        h(i, j) += (4 * m_T(i, j, k, l) + 8 * m_T(i, k, j, l)) * x[k] * x[l];
        return h;
    }

private:
    ElasticityTensor<Real, N> m_T;
};

#endif /* end of include guard: QUARTIC_FORM_HH */
