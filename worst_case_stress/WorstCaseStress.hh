////////////////////////////////////////////////////////////////////////////////
// WorstCaseStress.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computation and differentiation of microstructure worst-case stress
//      measures (both local and global).
//
//      For now, ***only piecewise constant per-element quantities are used***
//      (strains and stresses must first be averaged over each element).
//
//      WARNING: Cbase is assumed constant throughout the structure (i.e. single
//      base material). Variable materials can be added, but it will require
//      more data storage (or complicate the code).
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/03/2015 21:28:18
////////////////////////////////////////////////////////////////////////////////
#ifndef WORSTCASESTRESS_HH
#define WORSTCASESTRESS_HH
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <LinearElasticity.hh>
#include <PeriodicHomogenization.hh>
#include <tuple>
#include <GlobalBenchmark.hh>

#include "../pattern_optimization/SDConversions.hh"

// Local alias of PeriodicHomogenization namespace.
namespace {
    namespace PH = PeriodicHomogenization;
}

template<size_t N>
using MinorSymmetricRank4Tensor = ElasticityTensor<Real, N, false>;

// Piecewise constant rank 4 tensor fields.
template<size_t N>
using MinorSymmetricRank4TensorField = std::vector<MinorSymmetricRank4Tensor<N>>;

// (Pointwise) Worst-case stress measure
// Constructed by worstCase*Stress() below.
// This supports three types of worst-case stress:
//      1) Frobenius norm
//      2) von Mises
//      3) Max stress
// (1) and (2) are identical except (2) measures deviatoric micro stress.
// This means that for (2) F maps macro stress to deviatoric micro stress, and
// Cbase is really Dev : Cbase where Dev is the linear stress -> deviatoric
// stress map (I_4 - 1/3 I_2 tr).
//
// The max stress case (3) is the same as (1), except rank 4 tensor
// eigenvectors are and their eigenvalue are found (instead of eigenstrains).
// The eigenvectors n are converted into the corresponding worst-case macro
// stresses n n^T, and everything else (objective evaluation, shape
// derivatives) is identical.
template<size_t N>
struct WorstCaseStress {
    using SMF = SymmetricMatrixField<Real, N>;

    // d(wc stress)/d(kl^th cell problem strain) on element e
    SymmetricMatrixValue<Real, N> sensitvityToCellStrain(size_t kl, size_t e) const {
        assert(kl < flatLen(N));
        auto result = Cbase.doubleContract(F.at(e).doubleContract(wcMacroStress(e)));
        result *= 2 * (Sh.doubleContract(wcMacroStress(e)))[kl];
        return result;
    }

    // d(wc stress)/d(kl^th cell problem strain) rank 2 tensor field
    void sensitvityToCellStrain(size_t kl, SMF &result) const {
        assert(kl < flatLen(N));
        result.resizeDomain(size());
        for (size_t e = 0; e < size(); ++e) {
            result(e)  = Cbase.doubleContract(F.at(e).doubleContract(wcMacroStress(e)));
            result(e) *= 2 * (Sh.doubleContract(wcMacroStress(e)))[kl];
        }
    }

    // Get microscopic stress measure on element i.
    // Note: this is actually the squared Frobenius/von Mises/max stress
    Real operator()(size_t i) const {
        auto micro = F.at(i).doubleContract(wcMacroStress(i));
        return micro.doubleContract(micro);
    };

    // Get microscopic stress measure field.
    // Note: this is actually the squared Frobenius/von Mises/max stress
    ScalarField<Real> stressMeasure() const {
        ScalarField<Real> result(wcMacroStress.domainSize());
        for (size_t i = 0; i < result.domainSize(); ++i)
            result[i] = (*this)(i);
        return result;
    }

    // Get Frobenius/von Mises/max stress
    ScalarField<Real> sqrtStressMeasure() const {
        ScalarField<Real> result(wcMacroStress.domainSize());
        for (size_t i = 0; i < result.domainSize(); ++i)
            result[i] = sqrt((*this)(i));
        return result;
    }

    // Note: micro deviatoric stress in von Mises case.
    SMF wcMicroStress() const {
        SMF result(wcMacroStress.domainSize());
        for (size_t i = 0; i < F.size(); ++i)
            result(i) = F[i].doubleContract(wcMacroStress(i));
        return result;
    }

    size_t size() const {
        assert(F.size() == wcMacroStress.domainSize());
        return F.size();
    }

    ////////////////////////////////////////////////////////////////////////////
    // Discrete shape derivatives (Lagrangian). Useful for forward-mode diff.
    ////////////////////////////////////////////////////////////////////////////
    template<class Sim>
    ScalarField<Real> deltaStressMeasure(const Sim &sim,
                const std::vector<typename Sim::VField> &w,
                const typename Sim::VField &delta_p) const {
        size_t numElements = size();
        assert(numElements == sim.mesh().numElements());

        auto deltaSh = PH::deltaHomogenizedComplianceTensor(sim, w, delta_p);
        auto delta_w = PH::deltaFluctuationDisplacements(   sim, w, delta_p);
        auto deltaG  = PH::deltaMacroStrainToMicroStrainTensors(sim, w, delta_w, delta_p);

        MinorSymmetricRank4TensorField<N> deltaF; deltaF.reserve(deltaG.size());
        for (size_t e = 0; e < numElements; ++e) {
            deltaF.push_back(Cbase.doubleContract(deltaG[e].doubleContract(Sh)));
            deltaF.back() += Cbase.doubleContract(G[e].doubleContract(deltaSh));
        }

        ScalarField<Real> delta_s(numElements);
        for (size_t e = 0; e < numElements; ++e) {
            delta_s[e] = 2 * F[e].doubleContract(wcMacroStress(e)).doubleContract(
                    deltaF[e].doubleContract(wcMacroStress(e)));
        }

        return delta_s;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Public data members
    ////////////////////////////////////////////////////////////////////////////
    ElasticityTensor<Real, N> Cbase, Sh;
    // Per-element macro->micro {stress, deviatoric stress} map.
    // See Worst Case Microstructure writeup for details.
    MinorSymmetricRank4TensorField<N> F, G;
    SMF wcMacroStress;
    std::vector<int>  eigAlgebraicMult;
    std::vector<Real> eigPrincipal, eigSecondary;
};

template<size_t N>
WorstCaseStress<N> worstCaseFrobeniusStress(
        ElasticityTensor<Real, N> Cbase,
        ElasticityTensor<Real, N> Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain)
{
    WorstCaseStress<N> result;
    result.Cbase = Cbase;
    result.Sh = Sh;
    result.G = m2mStrain;

    size_t numElems = m2mStrain.size();

    // macro -> micro stress tensor map
    result.F.reserve(numElems);
    for (size_t i = 0; i < numElems; ++i)
        result.F.emplace_back(Cbase.doubleContract(m2mStrain[i].doubleContract(Sh)));

    // worst-case macro stress
    result.wcMacroStress.resizeDomain(numElems);
    result.eigAlgebraicMult.resize(numElems);
    result.eigPrincipal.resize(numElems);
    result.eigSecondary.resize(numElems);
    for (size_t i = 0; i < numElems; ++i) {
        // Actually major symmetric, but conversion is not yet supported
        MinorSymmetricRank4Tensor<N> T = result.F[i].transpose().doubleContract(result.F[i]);
        SymmetricMatrixValue<Real, N> sigma;
        std::tie(sigma, result.eigPrincipal[i],
                 result.eigAlgebraicMult[i],
                 result.eigSecondary[i]) = T.maxEigenstrainMultiplicity();
        result.wcMacroStress(i) = sigma;
    }
    return result;
}

template<size_t N>
WorstCaseStress<N> worstCaseVonMisesStress(
        ElasticityTensor<Real, N> Cbase,
        ElasticityTensor<Real, N> Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain)
{
    WorstCaseStress<N> result;
    // TODO
    throw std::runtime_error("Deviatoric stress extractor unimplemented");
    return worstCaseFrobeniusStress(/* DeviatoricExtractor.doubleContract( */ Cbase /* ) */,
                                    Sh, m2mStrain);
}

template<size_t N>
WorstCaseStress<N> worstCaseMaxStress(
        ElasticityTensor<Real, N> Cbase,
        ElasticityTensor<Real, N> Sh,
        const MinorSymmetricRank4TensorField<N> &m2mStrain)
{
    WorstCaseStress<N> result;
    result.Cbase = Cbase;
    result.Sh = Sh;
    result.G = m2mStrain;

    size_t numElems = m2mStrain.size();

    // macro -> micro stress tensor map
    result.F.reserve(numElems);
    for (size_t i = 0; i < numElems; ++i)
        result.F.emplace_back(Cbase.doubleContract(m2mStrain[i].doubleContract(Sh)));

    // worst-case macro stress
    result.wcMacroStress.resizeDomain(numElems);
    for (size_t e = 0; e < numElems; ++e) {
        // Actually major symmetric, but conversion is not yet supported
        MinorSymmetricRank4Tensor<N> T = result.F[e].transpose().doubleContract(result.F[e]);

        throw std::runtime_error("Rank 4 Tensor Eigenvectors Unimplemented");
        // TODO
        // SymmetricMatrixValue<Real, N> sigma;
        // Real lambda;
        VectorND<N> n(VectorND<N>::Zero());
        // tie(n, lambda) = T[e].maxEigenvector();

        // Worst-case stress is outer product of max eigenvalue's eigenvector.
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i; j < N; ++j)
                result.wcMacroStress(e)(i, j) = n[i] * n[j];
        }
    }
    return result;
}

// Global worst-case objective of the form
// int j(s, x) dV
// Currently, j (and s) are considered piecewise constant per-element.
template<size_t N, class Integrand>
struct IntegratedWorstCaseObjective {
    using SMF = typename WorstCaseStress<N>::SMF;

    IntegratedWorstCaseObjective() { }
    template<class Mesh> IntegratedWorstCaseObjective(const Mesh &m, const WorstCaseStress<N>  &wcs) { setPointwiseWCS(m, wcs); }
    template<class Mesh> IntegratedWorstCaseObjective(const Mesh &m,       WorstCaseStress<N> &&wcs) { setPointwiseWCS(m, std::move(wcs)); }

    template<class Mesh> void setPointwiseWCS(const Mesh &m, const WorstCaseStress<N>  &wcs) { wcStress = wcs;            integrand.init(m, wcStress); m_updateEvalCache(m); }
    template<class Mesh> void setPointwiseWCS(const Mesh &m,       WorstCaseStress<N> &&wcs) { wcStress = std::move(wcs); integrand.init(m, wcStress); m_updateEvalCache(m); }

    // Evaluate objective by integrating over m
    Real evaluate() const { return m_evalCache; }

    ScalarField<Real> integrandValues() const {
        const size_t numElems = wcStress.size();
        ScalarField<Real> result(numElems);
        for (size_t i = 0; i < numElems; ++i)
            result(i) = integrand.j(wcStress(i), i);
        return result;
    }

    // Compute the Eulerian derivative of the integrand at each element given
    // the per-element Eulerian derivative of fluctuation strains, dot_we.
    ScalarField<Real> integrandEulerianDerivative(
            const std::vector<SMF> &dot_we) const {
        size_t numElems = wcStress.size();
        ScalarField<Real> result(numElems);
        result.clear();
        SMF tau;
        for (size_t kl = 0; kl < dot_we.size(); ++kl) {
            assert(dot_we[kl].domainSize() == numElems);
            tau_kl(kl, tau);
            for (size_t i = 0; i < numElems; ++i) {
                // Sum is only over the "upper triangle" of integrands
                Real shearDoubler = (kl >= N) ? 2.0 : 1.0;
                result[i] += shearDoubler * (tau(i).doubleContract(dot_we[kl](i)));
            }
        }
        return result;
    }

    // Sensitivity of global objective integrand to cell problem fluctuation
    // strains (rank 2 tensor field tau^kl in the writeup).
    void tau_kl(size_t kl, SMF &result) const {
        wcStress.sensitvityToCellStrain(kl, result);
        for (size_t e = 0; e < wcStress.size(); ++e)
            result(e) *= integrand.j_prime(wcStress(e), e);
    }

    // Per-boundary-element interpolant to be integrated against normal shape
    // velocity to evaluate the shape derivative. This can be interpreted as
    // the objective's steepest ascent normal velocity field:
    //      j - strain(lambda^pq) : C : [strain(w^pq) + e^pq]
    //              + jShapeDependenceViaCh
    template<class Sim>
    using StrainEnergyBdryInterpolant = Interpolant<Real, Sim::K - 1, 2 * Sim::Strain::Deg>;
    template<class Sim>
    auto gradient(const Sim &sim,
                  const std::vector<VectorField<Real, N>> &w) const
    -> std::vector<StrainEnergyBdryInterpolant<Sim>>
    {
        BENCHMARK_START_TIMER("WCS shape derivative");
        using Strain = typename Sim::Strain;
        using BdryStrain = Interpolant<typename Sim::SMatrix, Sim::K - 1, Strain::Deg>;

        const auto &mesh = sim.mesh();

        std::vector<VectorField<Real, N>> lambda;
        m_solveAdjointCellProblems(sim, lambda);

        std::vector<StrainEnergyBdryInterpolant<Sim>> result(mesh.numBoundaryElements());
        for (auto e : mesh.elements()) {
            if (!e.isBoundary()) continue;
            const auto C = e->E();
            for (auto f : e.interfaces()) {
                auto be = mesh.boundaryElement(f.boundaryEntity().index());
                if (!be) continue;
                auto &r = result.at(be.index());
                if (be->isInternal) { r = 0; continue; }
                r = integrand.j(wcStress(e.index()), e.index());
                for (size_t pq = 0; pq < flatLen(N); ++pq) {
                    // Sum is only over the "upper triangle" of integrands
                    Real shearDoubler = (pq >= Sim::N) ? 2.0 : 1.0;

                    Strain strain_lambda, strain_w;
                    BdryStrain strain_lambda_bdry, strain_u_bdry;
                    e->strain(e, lambda[pq], strain_lambda);
                    e->strain(e,      w[pq], strain_w);

                    restrictInterpolant(e, be, strain_lambda, strain_lambda_bdry);
                    restrictInterpolant(e, be, strain_w,      strain_u_bdry     );
                    strain_u_bdry += Sim::SMatrix::CanonicalBasis(pq);
                    r -= Interpolation<Sim::K - 1, 2 * Strain::Deg>::interpolant(
                            [&] (const VectorND<Simplex::numVertices(Sim::K - 1)> &p) {
                                return ( C.doubleContract(strain_lambda_bdry(p))
                                          .doubleContract(strain_u_bdry(p)) ) * shearDoubler;
                            });
                }
            }
        }

        // Add in term accounting for j's shape dependence via C^H
        auto gamma_dCh = jShapeDependenceViaCh(sim, w);
        for (auto be : mesh.boundaryElements())
            result[be.index()] += gamma_dCh[be.index()];

        BENCHMARK_STOP_TIMER("WCS shape derivative");

        return result;
    }

    // Compute the pointwise (Eulerian) derivative of integrand j(wcs) due to a
    // normal shape velocity vn:
    //      tau^kl : strain(wdot^kl[vn]) + 2 * j' * sigma : F^T : C^Base : G : dS^H[vn] : sigma
    template<class Sim, class _NormalShapeVelocity>
    ScalarField<Real> directIntegrandDerivative(const Sim &sim,
                                                const std::vector<VectorField<Real, N>> &w,
                                                const std::vector<VectorField<Real, N>> &dot_w,
                                                const _NormalShapeVelocity &vn) const {
        const auto &mesh = sim.mesh();
        ScalarField<Real> result(mesh.numElements());
        result.clear();
        // Compute the tau term contribution to each element
        SMF tau, dot_w_kl_strain;
        for (size_t kl = 0; kl < dot_w.size(); ++kl) {
            tau_kl(kl, tau);
            // Sum is only over the "upper triangle" of integrands
            Real shearDoubler = (kl >= Sim::N) ? 2.0 : 1.0;
            dot_w_kl_strain = sim.averageStrainField(dot_w[kl]);
            for (auto e : mesh.elements()) {
                result[e.index()] += shearDoubler * (tau(e.index())
                              .doubleContract(dot_w_kl_strain(e.index())));
            }
        }

        // Compute variation of Ch
        auto sdCh = PH::homogenizedElasticityTensorGradient(w, sim);
        ElasticityTensor<Real, N> dCh_vn;
        for (auto be : mesh.boundaryElements()) {
            const auto &vnb = vn[be.index()]; const auto &sdb = sdCh[be.index()];
            using  SDInterp = typename std::decay<decltype(sdCh[0])>::type;
            using NSVInterp = typename std::decay<decltype(  vn[0])>::type;
            dCh_vn += Quadrature<SDInterp::K, SDInterp::Deg + NSVInterp::Deg>::
                integrate([&] (const EvalPt<SDInterp::K> &pt) {
                    return vnb(pt) * sdb(pt);
                }, be->volume());
        }

        // Compute variation of Sh
        auto dSh_vn = wcStress.Sh.doubleDoubleContract(dCh_vn);
        dSh_vn *= -1.0;

        // Compute the dSh term contribution to each element
        for (auto e : mesh.elements()) {
            // sigma : F^T : C^Base : G : dSh[vn] : sigma
            // = (F : sigma) : (C^Base : G : dSh[vn] : sigma)
            //   MicroStress              dMicroStress
            size_t ei = e.index();
            auto microStress = wcStress.F[ei].doubleContract(wcStress.wcMacroStress(ei));
            auto dMicroStress = wcStress.Cbase.doubleContract(
                            wcStress.G[ei].doubleContract(
                                dSh_vn.doubleContract(wcStress.wcMacroStress(ei))));
            Real dShTerm = 2 * integrand.j_prime(wcStress(ei), ei)
                             * microStress.doubleContract(dMicroStress);
            result[e.index()] += dShTerm;
        }

        return result;
    }

    // Direct (non-adjoint) evaluation of objective's shape derivative on
    // normal velocity field vn:
    //  int_vol tau^kl : strain(wdot^kl) dV + int_bdry vn * j dA
    template<class Sim, class _NormalShapeVelocity>
    Real directDerivative(const Sim &sim,
                          const std::vector<VectorField<Real, N>> &w,
                          const _NormalShapeVelocity &vn) const {
        std::vector<VectorField<Real, N>> dot_w;
        PH::fluctuationDisplacementShapeDerivatives(sim, w, vn, dot_w /* ,
                WCStressOptimization::Config::get().projectOutNormalStress */);
        assert(dot_w.size() == w.size());

        const auto &mesh = sim.mesh();

        Real advectionTerm = 0;
        // Assumes j is piecewise constant!
        // TODO: implement boundary element->element map in FEMMesh to simplify
        for (auto e : mesh.elements()) {
            if (!e.isBoundary()) continue;
            Real jval = integrand.j(wcStress(e.index()), e.index());
            for (auto f : e.interfaces()) {
                auto be = mesh.boundaryElement(f.boundaryEntity().index());
                if (!be) continue;
                if (be->isInternal) continue;
                advectionTerm += vn[be.index()].integrate(be->volume()) * jval;
            }
        }

#if 0
        Real dotKLTerm = 0;
        SMF tau, dot_w_kl_strain;
        for (size_t kl = 0; kl < dot_w.size(); ++kl) {
            tau_kl(kl, tau);
            // Sum is only over the "upper triangle" of integrands
            Real shearDoubler = (kl >= Sim::N) ? 2.0 : 1.0;
            dot_w_kl_strain = sim.averageStrainField(dot_w[kl]);
            for (auto e : mesh.elements()) {
                dotKLTerm += shearDoubler * e->volume() * (tau(e.index())
                                .doubleContract(dot_w_kl_strain(e.index())));
            }
        }

        // std::cout << dotKLTerm << '\t' << advectionTerm << '\t';

        Real jDirectTerm = 0;
        auto int_dj_domega = jShapeDependenceViaCh(sim, w);
        using  SDInterp = typename std::decay<decltype(int_dj_domega[0])>::type;
        using NSVInterp = typename std::decay<decltype(           vn[0])>::type;
        static_assert(SDInterp::K == NSVInterp::K, "Invalid boundary interpolant simplex dimension");

        for (auto be : mesh.boundaryElements()) {
            const auto &vnb = vn[be.index()];
            const auto &sdb = int_dj_domega[be.index()];
            jDirectTerm += Quadrature<SDInterp::K, SDInterp::Deg + NSVInterp::Deg>::
                integrate([&] (const VectorND<be.numVertices()> &pt) {
                    return vnb(pt) * sdb(pt);
                }, be->volume());
        }
        return advectionTerm + dotKLTerm + jDirectTerm;
#endif
        ScalarField<Real> dj = directIntegrandDerivative(sim, w, dot_w, vn);
        Real result = advectionTerm;
        for (auto e : mesh.elements())
            result += dj[e.index()] * e->volume();
        return result;
    }

    // Partial ``shape derivative'' of objective J, considering only the effect
    // of the changing homogenized elaticity tensor. In term of the writeup,
    // this is:
    //      gamma : dC^H[v]
    template<class Sim>
    auto jShapeDependenceViaCh(const Sim &sim,
                               const std::vector<VectorField<Real, N>> &w) const
    -> std::vector<StrainEnergyBdryInterpolant<Sim>>
    {
        auto gamma = dJ_dCH(sim);
        size_t numBE = sim.mesh().numBoundaryElements();
        std::vector<StrainEnergyBdryInterpolant<Sim>> result(numBE);

        auto gradCh = PH::homogenizedElasticityTensorGradient(w, sim);

        // Compute each nodal value of the interpolant.
        for (size_t i = 0; i < numBE; ++i) {
            auto &r = result[i];
            const auto &gCh = gradCh[i];
            using GChInterp = typename std::decay<decltype(gCh)>::type;
            using   RInterp = typename std::decay<decltype(  r)>::type;
            static_assert(RInterp::size() == GChInterp::size(),
                         "Interpolant node count mismatch");
            for (size_t n = 0; n < r.size(); ++n)
                r[n] = gamma.quadrupleContract(gCh[n]);
        }

        return result;
    }

    // Compute objective's partial derivative with respect to the homogenized
    // *elasticity* tensor C^H. This is rank-4 tensor "gamma" from the writeup:
    // dJ/dC^H =
    // 2 * int_omega (j') (G^T : C^base : F : sigma^*) otimes sigma^* dV :: dS^H[v]
    // := -D :: dS^H[v]                 (D = -2 * int_omega (j')...)
    // =   D :: (S^H : dC^H[v] : S^H)
    // = (S^H : D : S^H) :: dC^H[v]     (S^H is major symmetric)
    template<class Sim>
    MinorSymmetricRank4Tensor<N> dJ_dCH(const Sim &sim) const {
        MinorSymmetricRank4Tensor<N> D;
        for (size_t col = 0; col < flatLen(N); ++col) {
            auto c = D.DColAsSymMatrix(col);
            for (auto e : sim.mesh().elements()) {
                size_t ei = e.index();
                SymmetricMatrixValue<Real, N> contrib =
                    wcStress.G[ei].transpose().doubleContract(
                        wcStress.Cbase.doubleContract(
                            wcStress.F[ei].doubleContract(
                                wcStress.wcMacroStress(ei))));
                contrib *= e->volume() * wcStress.wcMacroStress(ei)[col] *
                           integrand.j_prime(wcStress(ei), ei);
                c += contrib;
            }
        }
        D *= -2.0;
        return wcStress.Sh.doubleDoubleContract(D);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Discrete shape derivative of the integrated worst-case stress objective
    // under mesh vertex perturbation delta_p (forward mode).
    // J is an integral of a piecewise-constant per-element field, so it's easy
    // to differentiate:
    //      J = int_omega j(s) dV = sum_e j(s)_e vol(e)
    //     dJ = sum_e [j(s)_e dvol(e) + j'(s)_e ds_e vol(e)]
    ////////////////////////////////////////////////////////////////////////////
    template<class Sim>
    Real deltaJ(const Sim &sim,
                const std::vector<typename Sim::VField> &w,
                const typename Sim::VField &delta_p) const {
        // Homogenized elasticity tensor change term: gamma :: dCh[v]
        auto deltaCh = PH::deltaHomogenizedElasticityTensor(sim, w, delta_p);
        Real dJ = dJ_dCH(sim).quadrupleContract(deltaCh);
        const auto &mesh = sim.mesh();

        // int tau_kl : delta strain(w^kl) dV
        auto delta_w = PH::deltaFluctuationDisplacements(sim, w, delta_p);
        {
            SMF tau;
            for (size_t kl = 0; kl < delta_w.size(); ++kl) {
                // Sum is only over the "upper triangle" of integrands
                Real shearDoubler = (kl >= N) ? 2.0 : 1.0;
                tau_kl(kl, tau);
                auto delta_we = sim.deltaAverageStrainField(w[kl], delta_w[kl], delta_p);
                for (auto e : mesh.elements())
                    dJ += tau(e.index()).doubleContract(delta_we(e.index())) * e->volume() * shearDoubler;
            }
        }

        // Volume dilation term: int j div v dV
        std::vector<VectorND<N>> cornerPerturbations;
        for (auto e : mesh.elements()) {
            size_t ei = e.index();
            sim.extractElementCornerValues(e, delta_p, cornerPerturbations);
            dJ += integrand.j(wcStress(ei), ei) * e->relativeDeltaVolume(cornerPerturbations) * e->volume();
        }

        return dJ;
    }

    // Compute the linear functional (one-form) dJ[v], represented as a
    // per-vertex vector field to be dotted with a **periodic** mesh vertex
    // velocity.
    //      dJ[v] = sum_i <delta_j[i], v[i]>
    //      (delta_j[i][c] = partial_derivative(J, v[i][c]))
    // where the sum is over all mesh vertices. (Vertices on the periodic
    // boundary get "half" contributions)
    //
    // -delta_j can be interpreted as a steepest descent direction in the
    // vertex position space R^(N*|v|), but it is not a good descent direction
    // for non-uniform meshes. For a better, mesh-independent direction, one
    // should use a different metric (E.g. ||v||_M = int <v(x), v(x)> dx = v[i]
    // M_ij v[j], where M is the deg 1 mass matrix. Or for Newton's method, the
    // Hessian.).
    template<class Sim>
    ScalarOneForm<Sim::N>
    adjointDeltaJ(const Sim &sim, const std::vector<typename Sim::VField> &w) const {
        // Dilation and delta strain terms
        const auto mesh = sim.mesh();
        size_t nv = mesh.numVertices();
        ScalarOneForm<Sim::N> delta_j(nv);
        delta_j.clear();

        BENCHMARK_START_TIMER_SECTION("Adjoint Cell Problem");
        BENCHMARK_START_TIMER("Compute Tau_kl");
        // Cache tau_pq and adjoint solutions
        std::vector<SMF> tau(flatLen(N));
        for (size_t pq = 0; pq < flatLen(N); ++pq)
            tau_kl(pq, tau[pq]);
        BENCHMARK_STOP_TIMER("Compute Tau_kl");
        std::vector<VectorField<Real, N>> lambda;
        m_solveAdjointCellProblems(sim, tau, lambda);
        BENCHMARK_STOP_TIMER_SECTION("Adjoint Cell Problem");

        BENCHMARK_START_TIMER_SECTION("Dilation and Delta strain");
        for (auto e : mesh.elements()) {
            // j contribution to dilation integrand
            Real dilationIntegrand = integrand.j(wcStress(e.index()), e.index());
            const auto C = e->E();
            for (size_t pq = 0; pq < flatLen(N); ++pq) {
                // Sum is only over the "upper triangle" of integrands
                Real shearDoubler = (pq >= Sim::N) ? 2.0 : 1.0;
                using Strain = typename Sim::Strain;
                Strain strain_lambda, strain_u;
                e->strain(e, lambda[pq], strain_lambda);
                e->strain(e,      w[pq], strain_u);
                strain_u += Sim::SMatrix::CanonicalBasis(pq);

                using Stress = typename Sim::Stress;
                Stress stress_lambda, stress_u;
                for (size_t i = 0; i < Stress::numNodalValues; ++i) {
                    stress_lambda[i] = C.doubleContract(strain_lambda[i]);
                    stress_u     [i] = C.doubleContract(strain_u     [i]);
                }

                // pq^th contribution to dilation integrand.
                dilationIntegrand -= Quadrature<N, 2 * Strain::Deg>::integrate(
                        [&] (const EvalPt<N> &p) {
                    return strain_lambda(p).doubleContract(stress_u(p));
                }) * shearDoubler;

                // Delta strain term
                for (auto n : e.nodes()) {
                    auto gradPhi_n = e->gradPhi(n.localIndex());

                    Interpolant<VectorND<N>, N, Stress::Deg> glam_functional;
                    glam_functional = tau[pq](e.index()).contract(w[pq](n.index()));
                    for (size_t i = 0; i < glam_functional.size(); ++i) {
                        glam_functional[i] -= stress_lambda[i].contract(w[pq](n.index()))
                                            + stress_u[i].contract(lambda[pq](n.index()));
                    }

                    for (auto v_m : e.vertices()) {
                        auto gradLam_m = e->gradBarycentric().col(v_m.localIndex());
                        Interpolant<Real, N, Strain::Deg> mu_grad_lam_m;
                        for (size_t i = 0; i < mu_grad_lam_m.size(); ++i)
                            mu_grad_lam_m[i] = gradLam_m.dot(glam_functional[i]);
                        delta_j(v_m.index()) -= Quadrature<N, 2 * Strain::Deg>::
                            integrate([&](const EvalPt<N> &pt) { return
                                    (mu_grad_lam_m(pt) * gradPhi_n(pt)).eval();
                                }, e->volume() * shearDoubler);
                    }
                }
            }
            for (auto v : e.vertices()) {
                delta_j(v.index()) += dilationIntegrand * e->volume()
                    * e->gradBarycentric().col(v.localIndex());
            }
        }
        BENCHMARK_STOP_TIMER_SECTION("Dilation and Delta strain");

        // Gamma term
        BENCHMARK_START_TIMER_SECTION("Gamma Term");
        {
            MinorSymmetricRank4Tensor<N> gamma = dJ_dCH(sim);
            auto dCh = PH::homogenizedElasticityTensorDiscreteDifferential(w, sim);
            ScalarOneForm<N> gtermNew = compose([&](const ElasticityTensor<Real, N> &e) { return e.quadrupleContract(gamma); }, dCh);
            delta_j += gtermNew;

            // using ETensorSD = PH::BEHTensorGradInterpolant<Sim>;
            // using GTermSD   = Interpolant<Real, ETensorSD::K, ETensorSD::Deg>;
            // std::vector<ETensorSD> sdCh = PH::homogenizedElasticityTensorGradient(w, sim);
            // std::vector<GTermSD> gammaTermFunctional(sdCh.size());
            // for (size_t i = 0; i < sdCh.size(); ++i) {
            //     for (size_t j = 0; j < ETensorSD::size(); ++j)
            //         gammaTermFunctional[i][j] = gamma.quadrupleContract(sdCh[i][j]);
            // }
            // auto gterm_bdry = SDConversions::diff_bdry_from_nsv_functional(gammaTermFunctional, mesh);
            // ScalarOneForm<N> gtermOld(mesh.numVertices());
            // gtermOld.clear();
            // for (auto bv : mesh.boundaryVertices())
            //     gtermOld(bv.volumeVertex().index()) += gterm_bdry(bv.index());
            // MSHFieldWriter writer("debug_gamma.msh", mesh);
            // writer.addField("gamma term new", gtermNew.asVectorField());
            // writer.addField("gamma term old", gtermOld.asVectorField());
        }
        BENCHMARK_STOP_TIMER_SECTION("Gamma Term");

        return delta_j;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Public data members
    ////////////////////////////////////////////////////////////////////////////
    WorstCaseStress<N> wcStress;
    Integrand integrand;

private:
    // TODO: change this to accept a BaseCellOperations class
    template<class Sim>
    void m_solveAdjointCellProblems(const Sim &sim,
                    std::vector<VectorField<Real, N>> &lambda_kl) const {
        size_t numCellProblems = flatLen(N);
        lambda_kl.clear(); lambda_kl.reserve(numCellProblems);
        SMF tau;
        for (size_t kl = 0; kl < numCellProblems; ++kl) {
            tau_kl(kl, tau);
            lambda_kl.push_back(sim.solve(sim.perElementStressFieldLoad(tau)));
        }
    }

    // Same as above, but with pre-computed tau_kl
    template<class Sim>
    void m_solveAdjointCellProblems(const Sim &sim,
                    const std::vector<SMF> &tau,
                    std::vector<VectorField<Real, N>> &lambda_kl) const {
        size_t numCellProblems = flatLen(N);
        lambda_kl.clear(); lambda_kl.reserve(numCellProblems);
        assert(tau.size() == numCellProblems);
        for (size_t kl = 0; kl < numCellProblems; ++kl)
            lambda_kl.push_back(sim.solve(sim.perElementStressFieldLoad(tau[kl])));
    }

    Real m_evalCache = 0.0;
    template<class Mesh> void m_updateEvalCache(const Mesh &m) {
        assert(m.numElements() == wcStress.size());
        m_evalCache = 0;
        for (auto e : m.elements())
            m_evalCache += integrand.j(wcStress(e.index()), e.index()) * e->volume();
    }
};

// int_omega worst_case_stress dV
// i.e. j(s, x) = s -> j' = 1.
struct WCStressIntegrandTotal {
    template<size_t N, class Mesh>
    void init(const Mesh &/* m */, const WorstCaseStress<N> &/* wcs */) { }

    // Derivative of global objective integrand wrt worst case stress.
    static Real j(Real wcStress, size_t /* x_i */) { return wcStress; }
    static Real j_prime(Real /* wcStress */, size_t /* x_i */) { return 1; }
};

// int_omega worst_case_stress^p dV
// i.e. j(s, x) = s^p -> j' = p s^(p - 1).
struct WCStressIntegrandLp {
    template<size_t N, class Mesh>
    void init(const Mesh &/* m */, const WorstCaseStress<N> &/* wcs */) { }

    WCStressIntegrandLp() { }

    Real j(Real wcStress, size_t /* x_i */) const { return pow(wcStress, p); }
    // Derivative of global objective integrand wrt worst case stress.
    Real j_prime(Real wcStress, size_t /* x_i */) const { return p * pow(wcStress, p - 1); }
    Real p = 1.0;
};

// Global max worst-case objective:
// max_{x in omega} s(x)
// We do approximate gradient computation by treating this as an integrated
// objective int_omega j(s, x) where
// j(s, x) = s * w(x)
// j'(s, x) = w(x)
// w(x) = 1.0 / V if worst case stress at s is at the max.
// V = total volume of all regions at max stress.
struct WCStressIntegrandLinf {
    template<size_t N, class Mesh>
    void init(const Mesh &m, const WorstCaseStress<N> &wcs) {
        assert(wcs.size() == m.numElements());
        assert(wcs.size() > 0);

        std::vector<size_t> maxStressElements(1, 0);
        Real maxStress = wcs(0);
        for (auto e : m.elements()) {
            Real val = wcs(e.index());
            if (maxStress < val) {
                maxStressElements.assign(1, e.index());
                maxStress = val;
            }
            else if (maxStress == val) {
                maxStressElements.push_back(e.index());
            }
        }

        Real volumeAtMaxStress = 0;
        for (size_t e : maxStressElements)
            volumeAtMaxStress += m.element(e)->volume();

        weightFunction.assign(m.numElements(), 0.0);
        for (size_t e : maxStressElements)
            weightFunction.at(e) = m.element(e)->volume() / volumeAtMaxStress;
    }

    Real j(Real s, size_t x_i) const { return s * weightFunction.at(x_i); }
    Real j_prime(Real /* s */, size_t x_i) const { return weightFunction.at(x_i); }

    std::vector<Real> weightFunction;
};

// Take the pth root of a WCS objective function.
// Based on chain rule, all derivatives are weighted by
// 1/p f^((1 / p) - 1), where f = SubObjective::eval();
template<class SubObjective>
struct PthRootObjective : public SubObjective {
using Base = SubObjective;
    // Forward all constructor args to SubObjective
    template<typename... Args>
    PthRootObjective(Args&&... args)
        : Base(std::forward<Args>(args)...) { }

    Real evaluate() const {
        if (p == 1) return Base::evaluate();
        return pow(Base::evaluate(), 1.0 / p);
    }

    template<class Sim, size_t N>
    auto gradient(const Sim &sim,
                  const std::vector<VectorField<Real, N>> &w) const
    -> std::vector<typename Base::template StrainEnergyBdryInterpolant<Sim>> {
        if (p == 1) return Base::gradient(sim, w);

        auto result = Base::gradient(sim, w);
        Real scale = m_gradientScale();
        for (auto &beInterp : result)
            beInterp *= scale;
        return result;
    }

    template<class Sim, class _NormalShapeVelocity, size_t N>
    Real directDerivative(const Sim &sim,
                          const std::vector<VectorField<Real, N>> &w,
                          const _NormalShapeVelocity &vn) const {
        Real pder = Base::directDerivative(sim, w, vn);
        if (p == 1) return pder;
        return m_gradientScale() * pder;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Discrete shape derivative under mesh vertex perturbation delta_p
    // (forward mode).
    ////////////////////////////////////////////////////////////////////////////
    template<class Sim>
    Real deltaJ(const Sim &sim,
                const std::vector<typename Sim::VField> &w,
                const typename Sim::VField &delta_p) const {
        Real pder = Base::deltaJ(sim, w, delta_p);
        if (p == 1) return pder;
        return m_gradientScale() * pder;
    }

    template<class Sim>
    ScalarOneForm<Sim::N>
    adjointDeltaJ(const Sim &sim,
            const std::vector<typename Sim::VField> &w) const {
        ScalarOneForm<Sim::N> pder = Base::adjointDeltaJ(sim, w); 
        if (p == 1) return pder;
        pder *= m_gradientScale();
        return pder;
    }

    Real p = 1.0;
private:
    Real m_gradientScale() const {
        return (1.0 / p) * pow(Base::evaluate(), 1.0 / p - 1.0);
    }
};


#endif /* end of include guard: WORSTCASESTRESS_HH */
