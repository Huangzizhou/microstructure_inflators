#include "WCSObjective.hh"
#include "WCStressOptimizationIterate.hh"
#include "../pattern_optimization/BoundaryPerturbationInflator.hh"
#include <iostream>

namespace WCStressOptimization {

template<class Sim, class WCSObjective = PthRootObjective<IntegratedWorstCaseObjective<Sim::N, WCStressIntegrandLp>>>
class BoundaryPerturbationIterate : public Iterate<Sim, WCSObjective, true /* _BypassParameterVelocity */> {
    using Base = Iterate<Sim, WCSObjective, true>;
    static constexpr size_t N = Sim::N;
    using ETensor = typename Sim::ETensor;
public:
    BoundaryPerturbationIterate(BoundaryPerturbationInflator<N> &inflator,
            size_t nParams, const double *params, Objective<N> &fullObjective)
        : Base(inflator, nParams, params, fullObjective), m_inflator(inflator)
    {
        assert(m_inflator.numParameters() == nParams);
    }

    // This high-dimensional iterate needs special adjoint-method based
    // gradient operations to be tractable.
    ScalarField<Real> gradientWCS_adjoint() const {
        auto sd = this->m_wcs_objective.gradient(*(this->m_sim), this->w_ij);
        if (Config::get().useVtxNormalPerturbationGradientVersion) return m_inflator.gradientFromShapeDerivativeVertexNormalVersion(sd);
        else                                                       return m_inflator.gradientFromShapeDerivative(sd);
    }
    ScalarField<Real> gradp_JS() const {
        auto sd = this->shapeDerivativeJS();
        if (Config::get().useVtxNormalPerturbationGradientVersion) return m_inflator.gradientFromShapeDerivativeVertexNormalVersion(sd);
        else                                                       return m_inflator.gradientFromShapeDerivative(sd);
    }

    ScalarField<Real> gradp_vol() const {
        auto sd = this->shapeDerivativeVolume();
        if (Config::get().useVtxNormalPerturbationGradientVersion) return m_inflator.gradientFromShapeDerivativeVertexNormalVersion(sd);
        else                                                       return m_inflator.gradientFromShapeDerivative(sd);
    }

    ScalarField<Real> gradp_JVol() const {
        auto sd = this->shapeDerivativeJVol();
        if (Config::get().useVtxNormalPerturbationGradientVersion) return m_inflator.gradientFromShapeDerivativeVertexNormalVersion(sd);
        else                                                       return m_inflator.gradientFromShapeDerivative(sd);
    }

protected:
    const BoundaryPerturbationInflator<N> &m_inflator;
};

}
