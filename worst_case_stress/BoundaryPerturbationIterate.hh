#include "WCStressOptimizationIterate.hh"
#include "../pattern_optimization/BoundaryPerturbationInflator.hh"
#include <iostream>

namespace WCStressOptimization {

template<class Sim, class Objective = PthRootObjective<IntegratedWorstCaseObjective<Sim::N, WCStressIntegrandLp>>>
class BoundaryPerturbationIterate : public Iterate<Sim, Objective, true /* _BypassParameterVelocity */> {
    using Base = Iterate<Sim, Objective, true>;
    static constexpr size_t N = Sim::N;
    using ETensor = typename Sim::ETensor;
public:
    BoundaryPerturbationIterate(BoundaryPerturbationInflator<N> &inflator,
            size_t nParams, const double *params, const ETensor &targetS)
        : Base(inflator, nParams, params, targetS), m_inflator(inflator)
    {
        assert(m_inflator.numParameters() == nParams);
    }

    void writeDescription(std::ostream &os) const {
        os << "moduli:\t";
        this->C.printOrthotropic(os);
        os << "anisotropy:\t" << this->C.anisotropy() << std::endl;

        os << "printable:\t" << this->m_printable << std::endl;

        os << "JS:\t"     << this->evaluateJS()    << std::endl;
        os << "WCS:\t"    << this->evaluateWCS()   << std::endl;
        if (this->m_targetVol)
            os << "JVol:\t"   << this->evaluateJVol()  << std::endl;
        os << "Volume:\t" << this->mesh().volume() << std::endl;

        auto gradP = gradientWCS_adjoint();
        os << "||grad_p WCS||:\t" << gradP.norm() << std::endl;;

        gradP = gradp_JS();
        os << "||grad_p JS||:\t" << gradP.norm() << std::endl;

        if (this->m_targetVol) {
            gradP = gradp_JVol();
            os << "||grad_p Jvol||:\t" << gradP.norm() << std::endl;
        }
    }

    // This high-dimensional iterate needs special adjoint-method based
    // gradient operations to be tractable.
    ScalarField<Real> gradientWCS_adjoint() const {
        auto sd = this->m_objective.gradient(*(this->m_sim), this->w_ij);
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
