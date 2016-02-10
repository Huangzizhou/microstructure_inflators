////////////////////////////////////////////////////////////////////////////////
// WCSObjective.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Configures the objective used for worst-case-stress optimization.
//      The objective is a weighted combination of several sub-objectives:
//          JS:   compliance tensor fitting
//          WCS:  worst-case stress measure (configured separately by
//                                           WCStressOptimizationConfig.hh)
//          JVol: volume target fitting
//          LReg: Laplacian of shape boundary (for regularization)
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  02/10/2016 12:57:34
////////////////////////////////////////////////////////////////////////////////
#ifndef WCSOBJECTIVE_HH
#define WCSOBJECTIVE_HH

#include <ElasticityTensor.hh>

namespace WCStressOptimization {

template<size_t N>
struct Objective {
    using ETensor = ElasticityTensor<Real, N>;
    using SField = ScalarField<Real> ;

    Objective(const ETensor _targetS, Real _jsWeight, Real _wcsWeight, Real
            _jvolWeight, Real _laplacianRegularizationWeight)
        : targetS(_targetS), jsWeight(_jsWeight), wcsWeight(_wcsWeight), jvolWeight(_jvolWeight),
          laplacianRegularizationWeight(_laplacianRegularizationWeight)
    { }

    // Normalized weights
    Real  normalizedWCSWeight() const { return wcsWeight / initialWCS(); }
    Real   normalizedJSWeight() const { return jsWeight / targetSNormSq(); }
    Real normalizedJVolWeight() const { return jvolWeight / targetVolSq(); }

    void setInitialWCS(Real wcs)          { m_hasInitialWCS   = true; m_initialWCS = wcs; }
    bool hasInitialWCS()            const { return m_hasInitialWCS; }

    void setTargetVolume(Real tvol)       { m_hasTargetVolume = true;  m_targetVol = tvol; }
    bool hasTargetVolume()          const { return m_hasTargetVolume; }
    Real targetVolume()             const { m_assertTVol(); return m_targetVol; }

    // Evaluate full objective from sub-objective values
    Real eval(Real JS, Real WCS, Real JVol = -1.0) const {
        if (JVol == -1.0) m_assertNoTVol();
        else              m_assertTVol();
        Real JFull =  JS * normalizedJSWeight() +
                      WCS * normalizedWCSWeight();
        if (hasTargetVolume()) JFull += JVol * normalizedJVolWeight();
        return JFull;
    }

    // Evaluate full objective gradient from sub-objective gradients
    SField evalGradient(const SField &JS_p, const SField &WCS_p, const SField &JVol_p = SField()) {
        if (JVol_p.domainSize() > 0) m_assertTVol();
        else                         m_assertNoTVol();
        SField grad_p = JS_p * normalizedJSWeight() +
                        WCS_p * normalizedWCSWeight();
        if (hasTargetVolume()) grad_p += JVol_p * normalizedJVolWeight();
        return grad_p;
    }

    // Quantities for normalization
    Real    initialWCS() const { m_assertInitWCS(); return m_initialWCS; }
    Real targetSNormSq() const { return targetS.quadrupleContract(targetS); }
    Real   targetVolSq() const { m_assertTVol(); return m_targetVol * m_targetVol; }
    ////////////////////////////////////////////////////////////////////////////
    // Public data members
    ////////////////////////////////////////////////////////////////////////////
    ETensor targetS;
    Real jsWeight, wcsWeight, jvolWeight, laplacianRegularizationWeight;
private:
    ////////////////////////////////////////////////////////////////////////////
    // Private data members
    ////////////////////////////////////////////////////////////////////////////
    bool m_hasInitialWCS = false;
    Real m_initialWCS = 0;
    bool m_hasTargetVolume = false;
    Real m_targetVol = 0.0;

    void m_assertTVol() const {
        if (!hasTargetVolume()) throw std::runtime_error("No target volume set.");
    }

    void m_assertNoTVol() const {
        if (!hasTargetVolume()) throw std::runtime_error("Expected a JVol (target volume is set).");
    }

    void m_assertInitWCS() const {
        if (!hasInitialWCS()) throw std::runtime_error("Initial WCS was not set. Make sure a common objective instance is used by all iterates.");
    }
};

}

#endif /* end of include guard: WCSOBJECTIVE_HH */
