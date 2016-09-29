#ifndef OBJECTIVETERMJS_HH
#define OBJECTIVETERMJS_HH

#include "../ObjectiveTerm.hh"
#include "../IterateFactory.hh"
#include "../SDConversions.hh"

#include <Fields.hh>
#include <OneForm.hh>
#include <ElasticityTensor.hh>
#include <LinearIndexer.hh>

#include <PeriodicHomogenization.hh>
#include <MSHFieldWriter.hh>

#include <stdexcept>
#include <iostream>

namespace PatternOptimization {
namespace ObjectiveTerms {

// NLLS: 1/2 ||S^H - S^*||_F^2
// The individual residual components are entries in the upper triangle of
// flattened tensor (S^H - S^*), weighted so that the squared residual vector
// norm computes the full rank 4 tensor Frobenius norm.
//
// For more accurate estimation, we estimate the residual and evaluate the
// objective in terms of it.
template<class _Sim>
struct TensorFit : NLLSObjectiveTerm<_Sim::N> {
    static constexpr size_t N = _Sim::N;
    using  OForm = ScalarOneForm<N>;
    using SField = ScalarField<Real>;
    using ETensor = ElasticityTensor<Real, N>;
    using VField = VectorField<Real, N>;
    template<class _Iterate>
    TensorFit(const ETensor &targetS, const _Iterate &it) : m_baseCellOps(it.baseCellOps()) {
        const auto S = it.complianceTensor();
        m_diffS = S - targetS;

        auto dChVol  = m_baseCellOps.homogenizedElasticityTensorDiscreteDifferential();
        auto dChBdry = m_baseCellOps.diff_bdry_from_diff_vol(dChVol);
        auto dShBdry = compose([&](const ETensor &e) { ETensor result = S.doubleDoubleContract(e); result *= -1.0; return result; }, dChBdry);
        this->m_differential = compose([&](const ETensor &e) { return m_diffS.quadrupleContract(e); }, dShBdry);

        {
            auto dJVol = compose([&](const ETensor &e) { return -m_diffS.quadrupleContract(S.doubleDoubleContract(e)); }, dChVol);
            MSHFieldWriter volWriter("volfields.msh", m_baseCellOps.mesh());
            volWriter.addField("dJVol", dJVol.asVectorField(), DomainType::PER_NODE);

            MSHBoundaryFieldWriter bdryWriter("bdryfields.msh", m_baseCellOps.mesh());
            const auto &sim = it.simulator();
            SField isInternal(sim.mesh().numBoundaryElements());
            isInternal.clear();
            for (auto be : sim.mesh().boundaryElements())
                isInternal(be.index()) = be->isInternal ? 1.0 : 0.0;
            bdryWriter.addField("isInternal", isInternal, DomainType::PER_ELEMENT);
        }


        // Test the iterate against 

        using LI = LinearIndexer<ETensor>;
        for (size_t i = 0; i < LI::size(); ++i) {
            m_component_differentials.push_back(
                    compose([&](const ETensor &e) { return LI::index(e, i); }, dShBdry));
        }
    }

    virtual Real evaluate() const override { return 0.5 * m_diffS.quadrupleContract(m_diffS); }

    bool ignoringShear() const { return m_ignoreShear; }
    void setIgnoreShear(bool ignore) { m_ignoreShear = ignore; }

    static constexpr size_t numResiduals() { return LinearIndexer<ETensor>::size(); }

	// The (ij, kl)th residual (kl >= ij) for the nonlinear least squares (a
    // single term of the Frobenius distance). The terms are weighted so
    // that the squared norm of the residual vector corresponds to the
    // Frobenius norm of the rank 4 tensor difference S - S^*.
    virtual SField residual() const override {
        SField result(numResiduals());
        for (size_t r = 0; r < numResiduals(); ++r) {
            size_t ij, kl;
            LinearIndexer<ETensor>::linearIndexTo2D(r, ij, kl);
            assert(kl >= ij);
            Real weight = 1.0;
            if (kl != ij) weight *= sqrt(2); // Account for lower triangle
            if (ignoringShear()) {
                if (ij >= N) weight = 0.0; // Zero out shear components
                if (kl >= N) weight = 0.0; // Zero out shear components
            }
            else {
                if (ij >= N) weight *= sqrt(2); // Left shear doubler
                if (kl >= N) weight *= sqrt(2); // Right shear doubler
            }
            result[r] = weight * m_diffS.D(ij, kl);
        }

        return result;
    }

    // Derivative of residual(ij, kl) wrt parameter p:
    // d/dp (S_ijkl - target_ijkl) = d/dp S_ijkl = <gradS_ijkl, vn_p>
    // The terms are weighted in accordance with the residual weighting above.
    virtual Eigen::MatrixXd jacobian(const std::vector<VField> &bdrySVels) const override {
        const size_t np = bdrySVels.size();
        Eigen::MatrixXd result(numResiduals(), np);
        for (size_t p = 0; p < np; ++p) {
            for (size_t r = 0; r < numResiduals(); ++r) {
                size_t ij, kl;
                LinearIndexer<ETensor>::linearIndexTo2D(r, ij, kl);
                Real weight = 1.0;
                if (kl != ij) weight *= sqrt(2); // Account for lower triangle
                if (ignoringShear()) {
                    if (ij >= N) weight = 0.0; // Zero out shear components
                    if (kl >= N) weight = 0.0; // Zero out shear components
                }
                else {
                    if (ij >= N) weight *= sqrt(2); // Left shear doubler
                    if (kl >= N) weight *= sqrt(2); // Right shear doubler
                }
                result(r, p) = weight * m_component_differentials.at(r)[bdrySVels[p]];
            }
        }
        return result;
    }

    virtual void writeFields(MSHFieldWriter &writer) const override {
        try {
            VField differential(m_baseCellOps.mesh().numVertices());
            differential.clear();
            for (auto bv : m_baseCellOps.mesh().boundaryVertices()) {
                differential(bv.volumeVertex().index()) =
                    this->m_differential.asVectorField()(bv.index());
            }
            writer.addField("JS Differential", differential, DomainType::PER_NODE);

            auto bdryVel = m_baseCellOps.descent_from_diff_bdry(this->m_differential);
            VField xferBdryVel(m_baseCellOps.mesh().numVertices());
            xferBdryVel.clear();
            for (auto v : m_baseCellOps.mesh().vertices()) {
                auto bv = v.boundaryVertex();
                if (!bv) continue;
                xferBdryVel(v.index()) = bdryVel(bv.index());
            }
            writer.addField("JS Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);
        }
        catch (const std::exception &e) {
            std::cerr << "Couldn't write TensorFit fields due to exception: " << e.what() << std::endl;
        }
    }

    virtual ~TensorFit() { }

private:
    // Differentials (one-forms) of each component of the compliance tensor
    const BaseCellOperations<_Sim> &m_baseCellOps;
    std::vector<OForm> m_component_differentials;
    bool m_ignoreShear = false;
    ETensor m_diffS;
};

// Configuration to be applyed by iterate factory
template<class _Sim>
struct IFConfigTensorFit : public IFConfig {
    static constexpr size_t N = _Sim::N;
    template<class _Iterate>
    void configIterate(const std::unique_ptr<_Iterate> &it, ObjectiveTermNormalizations &normalizations) const {
        static_assert(_Iterate::_N == N, "Mismatch in problem dimensions.");
        auto tf = Future::make_unique<TensorFit<_Sim>>(targetS, *it);
        tf->setWeight(weight);

        // Default JS normalization is target tensor's squared Frobenius norm
        if (!normalizations.isSet("JS"))
            normalizations.set("JS", 1.0 / targetS.frobeniusNormSq());

        tf->setNormalization(normalizations["JS"]);
        tf->setIgnoreShear(ignoreShear);
        it->addObjectiveTerm("JS", std::move(tf));
    }
    Real weight = 1.0;
    ElasticityTensor<Real, N> targetS;
    bool ignoreShear = false;
};

}}

#endif /* end of include guard: OBJECTIVETERMJS_HH */
