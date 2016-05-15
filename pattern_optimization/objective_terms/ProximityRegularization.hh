#ifndef OBJECTIVETERMPROXIMITYREGULARIZATION_HH
#define OBJECTIVETERMPROXIMITYREGULARIZATION_HH

#include "../ObjectiveTerm.hh"

namespace PatternOptimization {
namespace ObjectiveTerms {

// alpha * (p - p0).(p - p0) / 2
template<size_t N>
struct ProximityRegularization : public NLLSObjectiveTerm<N> {
    using SField = ScalarField<Real>;
    using VField = VectorField<Real, N>;

    ProximityRegularization(const std::vector<Real> &p,
                            const std::vector<Real> &initialParams)
        : m_currParams(p), m_initParams(initialParams) { }

    virtual SField gradp(const std::vector<VField> &/*bdrySVels*/) const {
        SField result(m_currParams);
        result -= m_initParams;
        result *= this->m_weight;
        return result;
    }

    virtual Real evaluate() const { return (this->m_weight / 2.0) * (m_initParams.values() - m_currParams.values()).squaredNorm(); }

    virtual SField residual() const {
        size_t nParams = m_initParams.domainSize();
        SField result(nParams);
        Real sqrtWeight = std::sqrt(this->m_weight);
        for (size_t i = 0; i < nParams; ++i)
            result[i] = sqrtWeight * (m_currParams[i] - m_initParams[i]);
        return result;
    }

    virtual Eigen::MatrixXd jacobian(const std::vector<VField> &/*bdrySVels*/) const {
        size_t nParams = m_initParams.domainSize();
        Eigen::MatrixXd result(nParams, nParams);
        result.setIdentity();
        result *= std::sqrt(this->m_weight);
        return result;
    }

    virtual ~ProximityRegularization() { }

private:
    SField m_currParams, m_initParams;
};

// Configuration to be applied by iterate factory
struct IFConfigProximityRegularization : public IFConfig {
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        if (!normalizations.isSet("ProximityRegularization"))
            normalizations.set("ProximityRegularization", 1.0); // TODO? Base on parameter range?

        auto pr = Future::make_unique<ProximityRegularization<_Iterate::_N>>(it->params(), initParams);
        pr->setWeight(weight);

        pr->setNormalization(normalizations["ProximityRegularization"]);
        it->addObjectiveTerm("ProximityRegularization", std::move(pr));
    }
    Real weight = 0.0;
    std::vector<Real> initParams;
};


}}

#endif /* end of include guard: OBJECTIVETERMPROXIMITYREGULARIZATION_HH */
