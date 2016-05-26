#ifndef OBJECTIVETERMJS_HH
#define OBJECTIVETERMJS_HH

#include "../ObjectiveTerm.hh"
#include "../IterateFactory.hh"
#include "../SDConversions.hh"

#include <PeriodicHomogenization.hh>
#include <MSHFieldWriter.hh>

namespace PatternOptimization {
namespace ObjectiveTerms {

namespace PH = PeriodicHomogenization;

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
    using SField = ScalarField<Real>;
    using ETensor = ElasticityTensor<Real, N>;
    using VField = VectorField<Real, N>;
    template<class _Iterate>
    TensorFit(const ETensor &targetS, const _Iterate &it) : m_sim(it.simulator()) {
        const auto S = it.complianceTensor();
        const auto &w = it.fluctuationDisplacements();
        m_diffS = S - targetS;

        // Compute compliance tensor gradient from elasticity tensor:
        // dSh = -S : dCh : S
        using ETensorSD = PH::BEHTensorGradInterpolant<_Sim>;
        std::vector<ETensorSD> dSh = PH::homogenizedElasticityTensorGradient(w, it.simulator());
        for (auto &GS : dSh) {
            for (size_t n = 0; n < GS.size(); ++n)
                GS[n] = -S.doubleDoubleContract(GS[n]);
        }

        // Compute differential of objective, first as NSV linear functional
        using NSVFuncInterp = Interpolant<Real, ETensorSD::K, ETensorSD::Deg>;
        std::vector<NSVFuncInterp> dJS(dSh.size());
        for (size_t i = 0; i < dSh.size(); ++i) {
            const auto &GS = dSh[i];
            for (size_t n = 0; n < GS.size(); ++n)
                dJS[i][n] = m_diffS.quadrupleContract(GS[n]);
        }
        this->m_differential = SDConversions::diff_bdry_from_nsv_functional(dJS, it.mesh());

        // compute each compliance tensor component's differentials
        std::vector<NSVFuncInterp> dShComponent(dSh.size());
        for (size_t r = 0; r < numResiduals(); ++r) {
            size_t ij, kl;
            residual2DIndexFrom1D(r, ij, kl);
            for (size_t bei = 0; bei < dSh.size(); ++bei) {
                const auto &GS = dSh[bei];
                for (size_t n = 0; n < GS.size(); ++n)
                    dShComponent[bei][n] = GS[n].D(ij, kl);
            }
            m_component_differentials.push_back(
                    SDConversions::diff_bdry_from_nsv_functional(dShComponent, it.mesh()));
        }
    }

    virtual Real evaluate() const override { return 0.5 * m_diffS.quadrupleContract(m_diffS); }

    bool ignoringShear() const { return m_ignoreShear; }
    void setIgnoreShear(bool ignore) { m_ignoreShear = ignore; }

    static constexpr size_t numResiduals() {
        return (flatLen(N) * (flatLen(N) + 1)) / 2;
    }

    // Compute the flattened tensor row and column index corresponding to a
    // residual component. Inverse of the following map:
    //    r(ij, kl) = kl + flatlen * ij - (ij * (ij + 1)) / 2
    // Not a closed form inverse, but likely faster than sqrt version anyway
    static void residual2DIndexFrom1D(const size_t r, size_t &ij, size_t &kl) {
        assert(r < numResiduals());
        size_t ri = r;
        kl = flatLen(N) + 1; // invalid
        for (ij = 0; ij < flatLen(N); ++ij) {
            size_t rowSize = flatLen(N) - ij;
            if (ri < rowSize) { kl = ri + ij; break;}
            ri -= rowSize;
        }
        assert((ij < flatLen(N)) && (kl < flatLen(N)));
    }

	// The (ij, kl)th residual (kl >= ij) for the nonlinear least squares (a
    // single term of the Frobenius distance). The terms are weighted so
    // that the squared norm of the residual vector corresponds to the
    // Frobenius norm of the rank 4 tensor difference S - S^*.
    virtual SField residual() const override {
        SField result(numResiduals());
        for (size_t r = 0; r < numResiduals(); ++r) {
            size_t ij, kl;
            residual2DIndexFrom1D(r, ij, kl);
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
                residual2DIndexFrom1D(r, ij, kl);
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
                result(r, p) = weight * m_component_differentials.at(r).innerProduct(bdrySVels[p]);
            }
        }
        return result;
    }

    virtual void writeFields(MSHFieldWriter &writer) const override {
        auto bdryVel = SDConversions::descent_from_diff_bdry(this->m_differential, m_sim);
        VField xferBdryVel(m_sim.mesh().numVertices());
        xferBdryVel.clear();
        for (auto v : m_sim.mesh().vertices()) {
            auto bv = v.boundaryVertex();
            if (!bv) continue;
            xferBdryVel(v.index()) = bdryVel(bv.index());
        }
        writer.addField("JS Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);
    }

    virtual ~TensorFit() { }

private:
    // Differentials (one-forms) of each component of the compliance tensor
    const _Sim &m_sim;
    std::vector<VField> m_component_differentials;
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
