////////////////////////////////////////////////////////////////////////////////
// WCStressOptimizationIterate.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Encapsulates the state of a worst-case stress pattern optimization
//      iterate and provides objective/gradient/etc.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/06/2015 23:16:57
////////////////////////////////////////////////////////////////////////////////
#ifndef WCSTRESSOPTIMIZATIONITERATE_HH
#define WCSTRESSOPTIMIZATIONITERATE_HH

#include <functional>
    
#include <PeriodicHomogenization.hh>
#include "../pattern_optimization/PatternOptimizationIterate.hh"
#include "WorstCaseStress.hh"

namespace WCStressOptimization {

template<class Sim>
class Iterate : public PatternOptimization::Iterate<Sim> {
    using Base = PatternOptimization::Iterate<Sim>;
    static constexpr size_t N = Sim::N;
    using ETensor = typename Sim::ETensor;
public:
    Iterate(ConstrainedInflator<N> &inflator, size_t nParams, const double *params,
            const ETensor &targetS)
        : Base(inflator, nParams, params, targetS, true /* Always keep w_ij */ )
    {
        // Worst case stress currently assumes that the base material is
        // constant, so we can read it off a single element.
        m_objective.setPointwiseWCS(
            worstCaseFrobeniusStress(m_sim->mesh().element(0)->E(), Base::S,
                PeriodicHomogenization::macroStrainToMicroStrainTensors(w_ij, *m_sim)));
    }

    // Evaluate the global worst case stress objective
    Real evaluateWCS() const { return m_objective.evaluate(m_sim->mesh()); }

    // Direct differentation version of gradient
    ScalarField<Real> gradientWCS_direct() const {
        ScalarField<Real> result(Base::m_params.size());
        for (size_t p = 0; p < result.domainSize(); ++p)
            result[p] = m_objective.directDerivative(*m_sim, w_ij, m_vn_p[p]);
        return result;
    }

    // Adjoint method version of gradient: evaluates linear functional form of
    // shape derivative
    ScalarField<Real> gradientWCS_adjoint() const {
        ScalarField<Real> result(Base::m_params.size());

        auto shapeDerivative = m_objective.gradient(*m_sim, w_ij);

        using  SDInterp = typename std::decay<decltype(shapeDerivative[0])>::type;
        using NSVInterp = typename std::decay<decltype(      m_vn_p[0][0])>::type;
        static_assert(SDInterp::K == NSVInterp::K,
                "Invalid boundary interpolant simplex dimension");

        // Boundary interpolant...
        for (size_t p = 0; p < result.domainSize(); ++p) {
            result[p] = 0;
            for (auto be : m_sim->mesh().boundaryElements()) {
                const auto &vn = m_vn_p[p][be.index()];
                const auto &sd = shapeDerivative[be.index()];
                result[p] += Quadrature<SDInterp::K, SDInterp::Deg + NSVInterp::Deg>::
                    integrate([&] (const VectorND<be.numVertices()> &pt) {
                        return vn(pt) * sd(pt);
                    }, be->volume());
            }
        }

        return result;
    }

protected:
    IntegratedWorstCaseObjective<N, WCStressIntegrandTotal> m_objective;
    using Base::m_sim;
    using Base::w_ij;
    using Base::m_vn_p;
};

}

#endif /* end of include guard: WCSTRESSOPTIMIZATIONITERATE_HH */
