#ifndef WCSOBJECTIVETERM_HH
#define WCSOBJECTIVETERM_HH

#include <ObjectiveTerm.hh>
#include <MSHFieldWriter.hh>

#include <Fields.hh>
#include <OneForm.hh>
#include <stdexcept>

#include "WorstCaseStress.hh"

namespace PatternOptimization {
namespace ObjectiveTerms {

template<class _Sim, class _WCSObjectiveType>
struct WorstCaseStress : ObjectiveTerm<_Sim::N> {
    using   Base = ObjectiveTerm<_Sim::N>;
    using  OForm = ScalarOneForm<_Sim::N>;
    using VField = typename _Sim::VField;
    template<class _Iterate>
    WorstCaseStress(const _Iterate &it,
                    Real globalObjectivePNorm, Real globalObjectiveRoot) : m_baseCellOps(it.baseCellOps()) {
        // Configure objective
        m_wcs_objective.integrand.p = globalObjectivePNorm;
        m_wcs_objective.p           = globalObjectiveRoot;

        const auto &mesh = it.simulator().mesh();

        // Worst case stress currently assumes that the base material is
        // constant, so we can read it off a single element.
        m_wcs_objective.setPointwiseWCS(mesh,
            worstCaseFrobeniusStress(mesh.element(0)->E(), it.complianceTensor(),
                m_baseCellOps.macroStrainToMicroStrainTensors()));

        // Compute and store WCS's boundary differential one-form
        // TODO: scale by numReflectedCells^(1.0 / p)
        m_diff_vol = m_wcs_objective.adjointDeltaJ(m_baseCellOps);
        this->m_differential = m_baseCellOps.diff_bdry_from_diff_vol(m_diff_vol);
    }

    // TODO: scale by numReflectedCells^(1.0 / p)
    virtual Real evaluate() const override { return m_wcs_objective.evaluate(); }

    virtual void writeFields(MSHFieldWriter &writer) const override {
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
        // writer.addField("Eigenvalue multiplicity", eigMult, DomainType::PER_ELEMENT);
        writer.addField("Eigenvalue relative distance", dist, DomainType::PER_ELEMENT);

        try {
            writer.addField("WCS Steepest Descent VVel",
                            m_baseCellOps.descent_from_diff_vol(m_diff_vol),
                            DomainType::PER_NODE);
        }
        catch (const std::exception &e) {
            std::cerr << "Couldn't write volume descent velocity: " << e.what() << std::endl;
        }

        auto bdryVel = m_baseCellOps.descent_from_diff_bdry(this->m_differential);
        VField xferBdryVel(m_baseCellOps.mesh().numVertices());
        xferBdryVel.clear();
        for (auto v : m_baseCellOps.mesh().vertices()) {
            auto bv = v.boundaryVertex();
            if (!bv) continue;
            xferBdryVel(v.index()) = bdryVel(bv.index());
        }
        writer.addField("WCS Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);
    }

    virtual void writeDescription(std::ostream &os, const std::string &name) const override {
        os << "Max Ptwise WCS:\t" << sqrt(m_wcs_objective.wcStress.stressMeasure().maxMag()) << std::endl;
        Base::writeDescription(os, name);
    }

    virtual ~WorstCaseStress() { }
private:
    OForm m_diff_vol; // per-volume-vertex differential
    const BaseCellOperations<_Sim> &m_baseCellOps;
    _WCSObjectiveType m_wcs_objective;
};

// Configuration to be applyed by iterate factory
template<class _Sim, class _WCSObjectiveType = PthRootObjective<IntegratedWorstCaseObjective<_Sim::N, WCStressIntegrandLp>>>
struct IFConfigWorstCaseStress : public IFConfig {
    static constexpr size_t N = _Sim::N;
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
        auto wcs = Future::make_unique<WorstCaseStress<_Sim, _WCSObjectiveType>>(
                *it, globalObjectivePNorm, globalObjectiveRoot);
        wcs->setWeight(weight);

        // WCS normalization is the initial worst case stress
        if (!normalizations.isSet("WCS"))
            normalizations.set("WCS", 1.0 / wcs->evaluate());

        wcs->setNormalization(normalizations["WCS"]);
        it->addObjectiveTerm("WCS", std::move(wcs));
    }

    Real weight = 1.0;
    Real globalObjectivePNorm = 1.0;
    Real globalObjectiveRoot = 1.0;
};

}}

#endif /* end of include guard: WCSOBJECTIVETERM_HH */
