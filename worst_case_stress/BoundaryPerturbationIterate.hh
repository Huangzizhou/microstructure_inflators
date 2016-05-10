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

    void writeDescription(std::ostream &os) const {
        writeIterateDescription(os, *this);
    }

    // This high-dimensional iterate needs special adjoint-method based
    // gradient operations to be tractable.
    ScalarField<Real> gradientWCS_adjoint() const {
        auto sd = this->m_wcs_objective.gradient(*(this->m_sim), this->w_ij);
        if (Config::get().useVtxNormalPerturbationGradientVersion) return m_inflator.gradientFromShapeDerivativeVertexNormalVersion(sd);
        else                                                       return m_inflator.gradientFromShapeDerivative(sd);
    }

    ScalarField<Real> gradientWCS_discrete_adjoint() const {
        typename Sim::VField sdb = this->steepestDescentBoundaryVelocity();
        sdb *= -1.0;
        assert(sdb.domainSize() == this->mesh().numBoundaryVertices());
        // extractParamsFromBoundaryValues needs a std::vector<VectorND>
        std::vector<VectorND<N>> bdryField;
        bdryField.reserve(sdb.domainSize());
        for (size_t i = 0; i < sdb.domainSize(); ++i) 
            bdryField.emplace_back(sdb(i));
        return m_inflator.extractParamsFromBoundaryValues(bdryField);
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

    ScalarField<Real> gradp_JFull() const {
        if (m_fullObjective.hasTargetVolume()) return m_fullObjective.evalGradient(gradp_JS(), gradientWCS_discrete_adjoint(), gradp_JVol());
        else                                   return m_fullObjective.evalGradient(gradp_JS(), gradientWCS_discrete_adjoint());
        // if (m_fullObjective.hasTargetVolume()) return m_fullObjective.evalGradient(gradp_JS(), gradientWCS_adjoint(), gradp_JVol());
        // else                                   return m_fullObjective.evalGradient(gradp_JS(), gradientWCS_adjoint());
    }

protected:
    const BoundaryPerturbationInflator<N> &m_inflator;
    using Base::m_fullObjective;
};

}
