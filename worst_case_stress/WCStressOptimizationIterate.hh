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
#include <MSHFieldWriter.hh>

#include "../pattern_optimization/PatternOptimizationIterate.hh"
#include "WorstCaseStress.hh"
#include "WCSObjective.hh"

namespace WCStressOptimization {

// Implemented outside a member method so that a subclasses' static overrides
// of the evaluation/gradient computations will be called.
// (E.g. BoundaryPerturbationIterate must call its own gradp_JS since it
// overrides the default parameter shape velocity gradient approach).
template<class _WCSIterate>
void writeIterateDescription(std::ostream &os, _WCSIterate &it) {
    os << "moduli:\t";
    it.elasticityTensor().printOrthotropic(os);
    os << "anisotropy:\t" << it.elasticityTensor().anisotropy() << std::endl;

    os << "printable:\t" << it.printable() << std::endl;

    BENCHMARK_START_TIMER("Evaluate JS/WCS");
    Real JS = it.evaluateJS(), WCS = it.evaluateWCS();
    os << "JS:\t"     << JS  << std::endl;
    os << "WCS:\t"    << WCS << std::endl;
    if (it.fullObjective().hasTargetVolume())
        os << "JVol:\t"   << it.evaluateJVol()  << std::endl;
    os << "Volume:\t" << it.mesh().volume() << std::endl;
    os << "Max Ptwise WCS:\t" << sqrt(it.wcsObjective().wcStress.stressMeasure().maxMag()) << std::endl;
    BENCHMARK_STOP_TIMER("Evaluate JS/WCS");

    auto   JS_p = it.gradp_JS();
    auto  WCS_p = it.gradientWCS_adjoint();
    ScalarField<Real> JVol_p(WCS_p.domainSize());
    JVol_p.clear();

    os << "||grad_p JS||:\t"  <<  JS_p.norm() << std::endl;
    os << "||grad_p WCS||:\t" << WCS_p.norm() << std::endl;

    if (it.fullObjective().hasTargetVolume()) {
        JVol_p = it.gradp_JVol();
        os << "||grad_p Jvol||:\t" << JVol_p.norm() << std::endl;
    }

    os << "Normalized WCS:\t" << WCS / it.fullObjective().initialWCS() << std::endl;
    os << "Normalized JS:\t" << it.evaluateJS() / it.fullObjective().targetSNormSq() << std::endl;
    if (it.fullObjective().hasTargetVolume())
        os << "Normalized JVol:\t" << it.evaluateJVol() / it.fullObjective().targetVolSq() << std::endl;

    // Output composite iterate stats.
    if (it.fullObjective().hasTargetVolume()) {
        os << "J_full:\t" << it.fullObjective().eval(JS, WCS, it.evaluateJVol()) << std::endl;
        os << "||grad_p J_full||:\t" << it.fullObjective().evalGradient(JS_p, WCS_p, JVol_p).norm() << std::endl;
    }
    else {
        os << "J_full:\t" << it.fullObjective().eval(JS, WCS) << std::endl;
        os << "||grad_p J_full||:\t" << it.fullObjective().evalGradient(JS_p, WCS_p).norm() << std::endl;
    }
    os << std::endl;
}

// _BypassParameterVelocity: hack to avoid computing the parameter velocity when
// there are too many parameters. This is only intended to be used by subclasses 
// that know how to compute objective gradients without explicitly forming the
// parameter velocity vectors (e.g. BoundaryPerturbationIterate)
template<class Sim, class WCSObjective = PthRootObjective<IntegratedWorstCaseObjective<Sim::N, WCStressIntegrandLp>>, bool _BypassParameterVelocity = false>
class Iterate : public PatternOptimization::Iterate<Sim, _BypassParameterVelocity> {
    using Base = PatternOptimization::Iterate<Sim, _BypassParameterVelocity>;
public:
    static constexpr size_t N = Sim::N;
    using ETensor = typename Sim::ETensor;
    using SField = ScalarField<Real>;

template<class _Inflator>
    // WARNING: updates fullObjective with init WCS on first construction (hack)
    Iterate(_Inflator &inflator, size_t nParams, const double *params,
            Objective<N> &fullObjective)
        : Base(inflator, nParams, params, fullObjective.targetS, true /* Always keep w_ij */ ),
          m_fullObjective(fullObjective)
    {
        if (m_fullObjective.hasTargetVolume())
            this->setTargetVolume(m_fullObjective.targetVolume());

        // Worst case stress currently assumes that the base material is
        // constant, so we can read it off a single element.
        m_wcs_objective.setPointwiseWCS(m_sim->mesh(),
            worstCaseFrobeniusStress(m_sim->mesh().element(0)->E(), Base::S,
                PeriodicHomogenization::macroStrainToMicroStrainTensors(w_ij, *m_sim)));

        if (!m_fullObjective.hasInitialWCS()) {
            m_fullObjective.setInitialWCS(evaluateWCS());
            std::cout << "Initial WCS: " << m_fullObjective.initialWCS() << std::endl;
        }
    }

    // Evaluate the global worst case stress objective
    Real evaluateWCS() const { return m_wcs_objective.evaluate(); }

    Real evaluateJFull() const {
        if (m_fullObjective.hasTargetVolume()) return m_fullObjective.eval(Base::evaluateJS(), evaluateWCS(), Base::evaluateJVol());
        else                                   return m_fullObjective.eval(Base::evaluateJS(), evaluateWCS());
    }

    // WARNING: paste this into any subclass that overrides gradient evaluation.
    SField gradp_JFull() const {
        if (m_fullObjective.hasTargetVolume()) return m_fullObjective.evalGradient(Base::gradp_JS(), gradientWCS_adjoint(), Base::gradp_JVol());
        else                                   return m_fullObjective.evalGradient(Base::gradp_JS(), gradientWCS_adjoint());
    }

    // Direct differentation version of gradient
    ScalarField<Real> gradientWCS_direct() const {
        ScalarField<Real> result(Base::m_params.size());
        for (size_t p = 0; p < result.domainSize(); ++p)
            result[p] = m_wcs_objective.directDerivative(*m_sim, w_ij, m_vn_p[p]);
        return result;
    }

    // Direct differentation version of gradient
    Real gradientWCS_direct_component(size_t p) const {
        return m_wcs_objective.directDerivative(*m_sim, w_ij, m_vn_p.at(p));
    }

    // Adjoint method version of gradient: evaluates linear functional form of
    // shape derivative
    ScalarField<Real> gradientWCS_adjoint() const {
        ScalarField<Real> result(Base::m_params.size());

        auto shapeDerivative = m_wcs_objective.gradient(*m_sim, w_ij);

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

    // WARNING: paste this into any subclass that overrides gradient/objective evaluation.
    void writeDescription(std::ostream &os) const {
        writeIterateDescription(os, *this);
    }

    void writeMeshAndFields(const std::string &name) const {
        MSHFieldWriter writer(name, m_sim->mesh());
        // typename Sim::VField outField;
        // for (size_t kl = 0; kl < flatLen(N); ++kl) {
        //     // Subtract off average displacements so that fields are comparable
        //     // across meshes.
        //     outField = w_ij[kl];
        //     outField -= outField.mean();
        //     writer.addField("w " + std::to_string(kl), outField);
        //     writer.addField("we " + std::to_string(kl), m_sim->averageStrainField(w_ij[kl]));
        // }
        ScalarField<Real> j = m_wcs_objective.integrandValues();

        // for (size_t p = 0; p < Base::m_params.size(); ++p) {
        //     auto prefix = "p" + std::to_string(p) + " ";
        //     std::vector<VectorField<Real, N>> dot_w;
        //     PeriodicHomogenization::fluctuationDisplacementShapeDerivatives(
        //             *m_sim, w_ij, m_vn_p[p], dot_w);
        //     std::vector<typename WCSObjective::SMF> wedot(flatLen(N));
        //     for (size_t kl = 0; kl < flatLen(N); ++kl) {
        //         outField = dot_w[kl];
        //         outField -= outField.mean();
        //         writer.addField(prefix + "wdot " + std::to_string(kl), outField);

        //         wedot[kl] = m_sim->averageStrainField(dot_w[kl]);
        //         writer.addField(prefix + "wedot " + std::to_string(kl), wedot[kl]);
        //     }
        //     writer.addField(prefix + "jdot", m_wcs_objective.integrandEulerianDerivative(wedot));

        //     // Validate tau_kl by calculating objective integrand at the offset
        //     // fluctuation displacements.
        //     auto offset_w = w_ij;
        //     Real delta = 1e-6; // see how size of this affects result...
        //     for (size_t kl = 0; kl < dot_w.size(); ++kl) {
        //         auto step = dot_w[kl];
        //         step *= delta;
        //         offset_w[kl] += step;
        //     }

        //     WCSObjective offsetObj(m_sim->mesh(),
        //             worstCaseFrobeniusStress(m_sim->mesh().element(0)->E(), m_wcs_objective.wcStress.Sh,
        //                 PeriodicHomogenization::macroStrainToMicroStrainTensors(offset_w, *m_sim)));
        //     ScalarField<Real> offsetDiff = offsetObj.integrandValues() - j;
        //     offsetDiff *= 1.0 / delta;
        //     writer.addField(prefix + "fd jdot", offsetDiff);

        //     // Same thing but with new compliance tensor
        //     // Why does offsetC differ from old compliance tensor????
        //     // Answer: it's because we're using the stress-like/displacement version
        //     // instead of the energy-like version. Though these versions are
        //     // equivalent for the true fluctuation displacements, they have
        //     // different behaviors for other displacements, and different
        //     // derivative behavior.
        //     // auto origChEnForm = PeriodicHomogenization::homogenizedElasticityTensor(w_ij, *m_sim);
        //     auto offsetCh = PeriodicHomogenization::homogenizedElasticityTensorDisplacementForm(offset_w, *m_sim);
        //     // auto offsetChEnForm = PeriodicHomogenization::homogenizedElasticityTensor(offset_w, *m_sim);

        //     // auto tmp = origChEnForm - m_wcs_objective.wcStress.Sh.inverse();
        //     // std::cout << "en form err:\t" << tmp.quadrupleContract(tmp) << std::endl;
        //     // tmp = offsetCh - m_wcs_objective.wcStress.Sh.inverse();
        //     // std::cout << "diff:\t" << tmp.quadrupleContract(tmp) << std::endl;
        //     // tmp = offsetChEnForm - m_wcs_objective.wcStress.Sh.inverse();
        //     // std::cout << "en form diff:\t" << tmp.quadrupleContract(tmp) << std::endl;

        //     offsetObj.setPointwiseWCS(
        //             worstCaseFrobeniusStress(m_sim->mesh().element(0)->E(), offsetCh.inverse(),
        //                 PeriodicHomogenization::macroStrainToMicroStrainTensors(offset_w, *m_sim)));
        //     offsetDiff = offsetObj.integrandValues() - j;
        //     offsetDiff *= 1.0 / delta;
        //     writer.addField(prefix + "corrected fd jdot", offsetDiff);

        //     auto gradEh =
        //         PeriodicHomogenization::homogenizedElasticityTensorGradient(w_ij, *m_sim);
        //     ETensor dC;
        //     for (auto be : m_sim->mesh().boundaryElements()) {
        //         const auto &vn = m_vn_p[p][be.index()];
        //         const auto &grad = gradEh[be.index()];

        //         using NSVInterp = typename std::decay<decltype(vn)>::type;
        //         using  SDInterp = typename std::decay<decltype(grad)>::type;
        //         static_assert(SDInterp::K == NSVInterp::K,
        //                 "Invalid boundary interpolant simplex dimension");

        //         dC += Quadrature<SDInterp::K, SDInterp::Deg + NSVInterp::Deg>::
        //             integrate([&] (const VectorND<be.numVertices()> &pt) {
        //                 return vn(pt) * grad(pt);
        //             }, be->volume());
        //     }
        //     offsetCh = m_wcs_objective.wcStress.Sh.inverse() + dC * delta;
        //     offsetObj.setPointwiseWCS(
        //             worstCaseFrobeniusStress(m_sim->mesh().element(0)->E(), offsetCh.inverse(),
        //                 PeriodicHomogenization::macroStrainToMicroStrainTensors(offset_w, *m_sim)));
        //     offsetDiff = offsetObj.integrandValues() - j;
        //     offsetDiff *= 1.0 / delta;
        //     writer.addField(prefix + "gradC corrected fd jdot", offsetDiff);
        //     
        // }

        writer.addField("WC Macro Stress", m_wcs_objective.wcStress.wcMacroStress);
        writer.addField("WC Micro Stress", m_wcs_objective.wcStress.wcMicroStress());

        writer.addField("Pointwise WCS", m_wcs_objective.wcStress.sqrtStressMeasure());
        writer.addField("j", j);
    }

    const WCSObjective &wcsObjective() const { return m_wcs_objective; }
    const WCStressOptimization::Objective<N> &fullObjective() const { return m_fullObjective; }

protected:
    WCStressOptimization::Objective<N> &m_fullObjective;
    WCSObjective m_wcs_objective;
    using Base::m_sim;
    using Base::w_ij;
    using Base::m_vn_p;
};

}

#endif /* end of include guard: WCSTRESSOPTIMIZATIONITERATE_HH */
