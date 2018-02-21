//
// Created by Davi Colli Tozoni on 1/9/18.
//

#ifndef STRESSOBJECTIVETERM_H
#define STRESSOBJECTIVETERM_H

#include <ObjectiveTerm.hh>
#include <MSHFieldWriter.hh>

#include <Fields.hh>
#include <OneForm.hh>
#include <stdexcept>

#include "MicroscopicStress.hh"

namespace PatternOptimization {
    namespace ObjectiveTerms {

        template<class _Sim, class _StressObjectiveType>
        struct StressTerm : ObjectiveTerm<_Sim::N> {
            using    Base = ObjectiveTerm<_Sim::N>;
            using   OForm = ScalarOneForm<_Sim::N>;
            using  VField = typename _Sim::VField;
            using SMatrix = typename _Sim::SMatrix;
            template<class _Iterate>
            StressTerm(const _Iterate &it, Real globalObjectivePNorm, Real globalObjectiveRoot, const std::string &measure)
                    : m_nonPeriodicCellOps(it.nonPeriodicCellOps())
            {
                // Configure objective
                m_stress_objective.integrand.p        = globalObjectivePNorm;
                m_stress_objective.p                  = globalObjectiveRoot;

                const auto &mesh = it.simulator().mesh();

                // We currently assumes that the base material is
                // constant, so we can read it off a single element.
                auto CBase = mesh.element(0)->E();
                if (measure == "frobenius") {
                    m_stress_objective.setPointwiseStress(mesh, MicroscopicFrobeniusStress<CBase.Dim, CBase.MajorSymmetry, _Sim>(CBase, it.simulator().averageStressField(m_nonPeriodicCellOps.displacement())));
                }
                else throw std::runtime_error("Unknown stress measure: " + measure);

                // Compute and store stress's boundary differential one-form
                m_diff_vol = m_stress_objective.adjointDeltaJ(m_nonPeriodicCellOps);
                this->m_differential = m_nonPeriodicCellOps.diff_bdry_from_diff_vol(m_diff_vol);
            }

            virtual Real evaluate() const override { return m_stress_objective.evaluate(); }

            virtual void writeFields(MSHFieldWriter &writer) const override {
                ScalarField<Real> j = m_stress_objective.integrandValues();

                writer.addField("Pointwise Stress", m_stress_objective.microStress.sqrtStressMeasure());
                // writer.addField("j", j);


#if 0
                try {
                    writer.addField("Stress Steepest Descent VVel",
                                    m_nonPeriodicCellOps.descent_from_diff_vol(m_diff_vol),
                                    DomainType::PER_NODE);
                }
                catch (const std::exception &e) {
                    std::cerr << "Couldn't write volume descent velocity: " << e.what() << std::endl;
                }
#endif

                auto bdryVel = m_nonPeriodicCellOps.descent_from_diff_bdry(this->m_differential);
                VField xferBdryVel(m_nonPeriodicCellOps.mesh().numVertices());
                xferBdryVel.clear();
                for (auto v : m_nonPeriodicCellOps.mesh().vertices()) {
                    auto bv = v.boundaryVertex();
                    if (!bv) continue;
                    xferBdryVel(v.index()) = bdryVel(bv.index());
                }
                writer.addField("Stress Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);

                auto dJ_field = m_nonPeriodicCellOps.descent_from_diff_vol(m_diff_vol);
                writer.addField("m_diff_vol", dJ_field, DomainType::PER_NODE);
            }

            virtual void writeDescription(std::ostream &os, const std::string &name) const override {
                os << "Ptwise stress:\t" << sqrt(m_stress_objective.microStress.stressMeasure().maxMag()) << std::endl;
                Base::writeDescription(os, name);
            }

            virtual ~StressTerm() { }
        private:
            OForm m_diff_vol; // per-volume-vertex differential
            const NonPeriodicCellOperations<_Sim> &m_nonPeriodicCellOps;
            _StressObjectiveType m_stress_objective;
        };

        // Configuration to be applied by iterate factory
        template<class _Sim, class _StressObjectiveType = PthRootObjective<IntegratedMicroscopicStressObjective<_Sim::N, MicroscopicStressIntegrandLp<_Sim>, _Sim>>>
        struct IFConfigMicroscopicStress : public IFConfig {
            static constexpr size_t N = _Sim::N;
            template<class _Iterate>
            void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
                static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
                BENCHMARK_START_TIMER_SECTION("Stress Term");
                auto stress_term = Future::make_unique<StressTerm<_Sim, _StressObjectiveType>>(
                        *it, globalObjectivePNorm, globalObjectiveRoot, measure);
                stress_term->setWeight(weight);

                // Stress normalization is the initial worst case stress
                if (!normalizations.isSet("Stress"))
                    normalizations.set("Stress", 1.0 / stress_term->evaluate());

                stress_term->setNormalization(normalizations["Stress"]);
                it->addObjectiveTerm("Stress", std::move(stress_term));
                BENCHMARK_STOP_TIMER_SECTION("Stress Term");
            }

            Real weight = 1.0;
            Real globalObjectivePNorm = 1.0;
            Real globalObjectiveRoot = 1.0;
            std::string measure = std::string("frobenius");
        };
    }
}

#endif //STRESSOBJECTIVETERM_H

