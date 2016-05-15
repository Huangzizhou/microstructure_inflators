#ifndef WCSOBJECTIVETERM_HH
#define WCSOBJECTIVETERM_HH

#include <ObjectiveTerm.hh>
#include <MSHFieldWriter.hh>

#include "WorstCaseStress.hh"

namespace PatternOptimization {
namespace ObjectiveTerms {

namespace PH = PeriodicHomogenization;

template<class _Sim, class _WCSObjectiveType>
struct WorstCaseStress : ObjectiveTerm<_Sim::N> {
    using Base = ObjectiveTerm<_Sim::N>;
    using VField = typename _Sim::VField;
    template<class _Iterate>
    WorstCaseStress(const _Iterate &it, const _Sim &sim) : m_sim(sim) {
        // Worst case stress currently assumes that the base material is
        // constant, so we can read it off a single element.
        m_wcs_objective.setPointwiseWCS(m_sim.mesh(),
            worstCaseFrobeniusStress(m_sim.mesh().element(0)->E(), it.complianceTensor(),
                PH::macroStrainToMicroStrainTensors(it.fluctuationDisplacements(), m_sim)));

        // Compute and store WCS's boundary differential one-form
        m_diff_vol = m_wcs_objective.adjointDeltaJ(m_sim, it.fluctuationDisplacements());
        
        this->m_differential  = SDConversions::diff_bdry_from_diff_vol(m_diff_vol, m_sim);
    }

    virtual Real evaluate() const { return m_wcs_objective.evaluate(); }

    virtual void writeFields(MSHFieldWriter &writer) const {
        ScalarField<Real> j = m_wcs_objective.integrandValues();

        writer.addField("Pointwise WCS", m_wcs_objective.wcStress.sqrtStressMeasure());
        // writer.addField("j", j);

        ScalarField<Real> eigPrincipal(m_wcs_objective.wcStress.eigPrincipal),
                          eigSecondary(m_wcs_objective.wcStress.eigSecondary),
                          eigMult(m_wcs_objective.wcStress.size()),
                          dist(m_wcs_objective.wcStress.size());
        for (size_t i = 0; i < eigMult.domainSize(); ++i) {
            eigMult[i] = Real(m_wcs_objective.wcStress.eigAlgebraicMult.at(i));
            dist[i] = (eigPrincipal[i] - eigSecondary[i]) / eigPrincipal[i];
        }
        // writer.addField("Principal eigenvalue", eigPrincipal, DomainType::PER_ELEMENT);
        // writer.addField("Secondary eigenvalue", eigSecondary, DomainType::PER_ELEMENT);
        writer.addField("Eigenvalue relative distance", dist, DomainType::PER_ELEMENT);

        writer.addField("Steepest Descent VVel",
                        SDConversions::descent_from_diff_vol(m_diff_vol, m_sim),
                        DomainType::PER_NODE);

        auto bdryVel = SDConversions::descent_from_diff_bdry(this->m_differential, m_sim);
        VField xferBdryVel(m_sim.mesh().numVertices());
        xferBdryVel.clear();
        for (auto v : m_sim.mesh().vertices()) {
            auto bv = v.boundaryVertex();
            if (!bv) continue;
            xferBdryVel(v.index()) = bdryVel(bv.index());
        }
        writer.addField("Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const {
        os << "Max Ptwise WCS:\t" << sqrt(m_wcs_objective.wcStress.stressMeasure().maxMag()) << std::endl;
        Base::writeDescription(os, name);
    }

    virtual ~WorstCaseStress() { }
private:
    VField m_diff_vol; // per-volume-vertex differential
    const _Sim &m_sim;
    _WCSObjectiveType m_wcs_objective;
};

// Configuration to be applyed by iterate factory
template<class _Sim, class _WCSObjectiveType = PthRootObjective<IntegratedWorstCaseObjective<_Sim::N, WCStressIntegrandLp>>>
struct IFConfigWorstCaseStress : public IFConfig {
    static constexpr size_t N = _Sim::N;
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
        auto wcs = Future::make_unique<WorstCaseStress<_Sim, _WCSObjectiveType>>(*it, it->simulator());
        wcs->setWeight(weight);

        // WCS normalization is the initial worst case stress
        if (!normalizations.isSet("WCS"))
            normalizations.set("WCS", 1.0 / wcs->evaluate());

        wcs->setNormalization(normalizations["WCS"]);
        it->addObjectiveTerm("WCS", std::move(wcs));
    }
    Real weight = 1.0;
};

}}

#endif /* end of include guard: WCSOBJECTIVETERM_HH */
