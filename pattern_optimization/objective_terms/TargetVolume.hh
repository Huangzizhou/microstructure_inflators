#ifndef OBJECTIVETERMJVOL_HH
#define OBJECTIVETERMJVOL_HH

#include "../ObjectiveTerm.hh"
#include "../IterateFactory.hh"

namespace PatternOptimization {
namespace ObjectiveTerms {

template<size_t N>
struct TargetVolume : public NLLSObjectiveTerm<N> {
    TargetVolume(Real tvol) : m_targetVol(tvol) { }
    
    void setTarget(v) { m_targetVol = v; }
    virtual ~TargetVolume() { }
private:
    Real m_targetVol;
};

// Configuration to be applyed by iterate factory
struct IFConfigTargetVolume : public IFConfig {
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        auto tv = Future::make_unique<ObjectiveTermJVol<_Iterate::_N>>(targetVolume, *it);
        tv->setWeight(weight);

        if (!normalizations.isSet("JVol"))
            normalizations.set("JVol", 1.0); // TODO? What about target vol = 0?

        it->addObjectiveTerm("JVol", std::move(tv));
        it->setNormalization(normalizations["JVol"]);
    }
    Real weight = 0.0;
    Real targetVolume = 0.0;
};

}}

#endif /* end of include guard: OBJECTIVETERMJVOL_HH */
