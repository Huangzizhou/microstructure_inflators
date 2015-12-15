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

    ElasticityTensor<Real, N> Cbase, Sh;
    // Per-element macro->micro {stress, deviatoric stress} map.
    // See Worst Case Microstructure writeup for details.
    MinorSymmetricRank4TensorField<N> F, G;
    SMF wcMacroStress;
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
    for (size_t i = 0; i < numElems; ++i) {
        // Actually major symmetric, but conversion is not yet supported
        MinorSymmetricRank4Tensor<N> T = result.F[i].transpose().doubleContract(result.F[i]);
        SymmetricMatrixValue<Real, N> sigma;
        Real lambda;
        std::tie(sigma, lambda) = T.maxEigenstrain();
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
    IntegratedWorstCaseObjective(const WorstCaseStress<N>  &wcs) { setPointwiseWCS(wcs); }
    IntegratedWorstCaseObjective(      WorstCaseStress<N> &&wcs) { setPointwiseWCS(std::move(wcs)); }

    void setPointwiseWCS(const WorstCaseStress<N>  &wcs) { wcStress = wcs; }
    void setPointwiseWCS(      WorstCaseStress<N> &&wcs) { wcStress = std::move(wcs); }

    // Evaluate objective by integrating over m
    template<class Mesh>
    Real evaluate(const Mesh &m) const {
        assert(m.numElements() == wcStress.size());
        Real result = 0;
        for (auto e : m.elements())
            result += Integrand::j(wcStress(e.index()), e.index()) * e->volume();
        return result;
    }

    ScalarField<Real> integrand() const {
        const size_t numElems = wcStress.size();
        ScalarField<Real> result(numElems);
        for (size_t i = 0; i < numElems; ++i)
            result(i) = Integrand::j(wcStress(i), i);
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
            result(e) *= Integrand::j_prime(wcStress(e), e);
    }

    // Per-boundary-element interpolant to be integrated against normal shape
    // velocity to evaluate the shape derivative. This can be interpreted as
    // the objective's steepest ascent normal velocity field:
    //      j - strain(lambda^pq) : C : [strain(w^pq) + e^pq]
    //              + jDirectShapeDependence
    template<class Sim>
    using StrainEnergyBdryInterpolant = Interpolant<Real, Sim::K - 1, 2 * Sim::Strain::Deg>;
    template<class Sim>
    auto gradient(const Sim &sim,
                  const std::vector<VectorField<Real, N>> &w) const
    -> std::vector<StrainEnergyBdryInterpolant<Sim>>
    {
        using Strain = typename Sim::Strain;
        using BdryStrain = Interpolant<typename Sim::SMatrix, Sim::K - 1, Strain::Deg>;

        const auto &mesh = sim.mesh();

        std::vector<VectorField<Real, N>> lambda;
        m_solveAdjointCellProblems(sim, lambda);

        std::vector<StrainEnergyBdryInterpolant<Sim>> result(mesh.numBoundaryElements());
        for (auto e : mesh.elements()) {
            if (!e.isBoundary()) continue;
            const auto C = e->E();
            for (size_t fi = 0; fi < e.numNeighbors(); ++fi) {
                auto be = mesh.boundaryElement(e.interface(fi).boundaryEntity().index());
                if (!be) continue;
                auto &r = result.at(be.index());
                if (be->isPeriodic) { r = 0; continue; }
                r = Integrand::j(wcStress(e.index()), e.index());
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

        // Add in term accounting for j's direct shape dependence
        auto int_dj_domega = jDirectShapeDependence(sim, w);
        for (auto be : mesh.boundaryElements())
            result[be.index()] += int_dj_domega[be.index()];
        
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
        PH::fluctuationDisplacementShapeDerivatives(sim, w, vn, dot_w);
        assert(dot_w.size() == w.size());

        const auto &mesh = sim.mesh();

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

        Real advectionTerm = 0;
        // Assumes j is piecewise constant!
        // TODO: implement boundary element->element map in FEMMesh to simplify
        for (auto e : mesh.elements()) {
            if (!e.isBoundary()) continue;
            Real jval = Integrand::j(wcStress(e.index()), e.index());
            for (size_t fi = 0; fi < e.numNeighbors(); ++fi) {
                auto be = mesh.boundaryElement(e.interface(fi).boundaryEntity().index());
                if (!be) continue;
                if (be->isPeriodic) continue;
                advectionTerm += vn[be.index()].integrate(be->volume()) * jval;
            }
        }

        std::cout << dotKLTerm << '\t' << advectionTerm << '\t';

        Real jDirectTerm = 0;
        auto int_dj_domega = jDirectShapeDependence(sim, w);
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

        return dotKLTerm + advectionTerm + jDirectTerm;
    }

    // Effect of j's direct dependence on the microstructure's shape.
    // int_omega dj/domega dx =
    // 2 * int_omega (j') G^T : C^base : F : sigma^* otimes sigma^* dV :: dS^H[v]
    // := -D :: dS^H[v]                 (D = -2 * int_omega (j')...)
    // =   D :: (S^H : dC^H[v] : S^H)
    // = (S^H : D : S^H) :: dC^H[v]     (S^H is major symmetric)
    //
    // TODO: make boundary velocity functional class to wrap these
    // interpolants (make evaluating/manipulating them easier).
    template<class Sim>
    auto jDirectShapeDependence(const Sim &sim,
                                const std::vector<VectorField<Real, N>> &w) const
    -> std::vector<StrainEnergyBdryInterpolant<Sim>>
    {
        const auto &mesh = sim.mesh();
        MinorSymmetricRank4Tensor<N> D;
        for (size_t col = 0; col < flatLen(N); ++col) {
            auto c = D.DColAsSymMatrix(col);
            for (auto e : mesh.elements()) {
                size_t ei = e.index();
                SymmetricMatrixValue<Real, N> contrib = 
                    wcStress.G[ei].transpose().doubleContract(
                        wcStress.Cbase.doubleContract(
                            wcStress.F[ei].doubleContract(
                                wcStress.wcMacroStress(ei))));
                contrib *= e->volume() * wcStress.wcMacroStress(ei)[col] *
                           Integrand::j_prime(wcStress(ei), ei);
                c += contrib;
            }
        }
        D *= -2.0;
        auto ShDSh = wcStress.Sh.doubleDoubleContract(D);

        size_t numBE = mesh.numBoundaryElements();
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
                r[n] = ShDSh.quadrupleContract(gCh[n]);
        }

        return result;
    }

    WorstCaseStress<N> wcStress;

private:
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
};

// int_omega worst_case_stress dV
// i.e. j(s, x) = s -> j' = 1.
struct WCStressIntegrandTotal {
    // Derivative of global objective integrand wrt worst case stress.
    static Real j(Real wcStress, size_t /* x_i */) { return wcStress; }
    static Real j_prime(Real /* wcStress */, size_t /* x_i */) { return 1; }
};

// int_omega worst_case_stress^p dV
// i.e. j(s, x) = s^p -> j' = p s^(p - 1).
template<int p>
struct WCStressIntegrandLp {
    static Real j(Real wcStress, size_t /* x_i */) { return pow(wcStress, p); }
    // Derivative of global objective integrand wrt worst case stress.
    static Real j_prime(Real wcStress, size_t /* x_i */) { return p * pow(wcStress, p - 1); }
};

#endif /* end of include guard: WORSTCASESTRESS_HH */
