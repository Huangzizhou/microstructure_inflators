////////////////////////////////////////////////////////////////////////////////
// MicroscopicStress.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Computation and differentiation of stress measures (both local and global).
//
//      For now, ***only piecewise constant per-element quantities are used***
//      (strains and stresses must first be averaged over each element).
//      Based on WorstCaseStress.hh
//
*/
//  Author:  Davi Colli Tozoni (dctozoni), davi.tozoni@nyu.edu
//  Company:  New York University
//  Created:  1/9/2018
////////////////////////////////////////////////////////////////////////////////

#ifndef MICROSCOPICSTRESS_H
#define MICROSCOPICSTRESS_H

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <LinearElasticity.hh>
#include <PeriodicHomogenization.hh>
#include <tuple>
#include <GlobalBenchmark.hh>
#include <VonMises.hh>
#include <Parallelism.hh>

#include "../pattern_optimization/SDConversions.hh"
#include "NonPeriodicCellOperations.hh"

template<size_t N>
using MinorSymmetricRank4Tensor = ElasticityTensor<Real, N, false>;

// Piecewise constant rank 4 tensor fields.
template<size_t N>
using MinorSymmetricRank4TensorField = std::vector<MinorSymmetricRank4Tensor<N>>;

// (Pointwise) Stress measure constructed by Microscopic*Stress() below.
// This supports three types of worst-case stress:
//      1) Frobenius norm (only one so far)
//      2) von Mises (TODO:implement)
//      3) Max stress (TODO:implement)
//      4) Stress trace (TODO:implement)
template<size_t N>
struct MicroscopicStress {
    using SMF = SymmetricMatrixField<Real, N>;

    // d(micro stress)/d(cell problem strain) on element e
    SymmetricMatrixValue<Real, N> sensitivityToCellStrain(size_t e) const {
        auto result = stressField[e].doubleContract(Cbase);
        result *= 2;
        return result;
    }

    // d(micro stress)/d(cell problem strain) rank 2 tensor field
    void sensitivityToCellStrain(SMF &result) const {
        result.resizeDomain(size());
        for (size_t e = 0; e < size(); ++e) {
            result(e)  = stressField[e].doubleContract(Cbase);
            result(e) *= 2;
        }
    }

    // Get microscopic stress measure on element i.
    // Note: this is actually the squared Frobenius stress
    Real operator()(size_t i) const {
        auto micro = stressField.at(i);
        return micro.doubleContract(micro);
    };

    // Get microscopic stress measure field.
    // Note: this is actually the squared Frobenius/von Mises/max stress
    ScalarField<Real> stressMeasure() const {
        ScalarField<Real> result(stressField.domainSize());
        for (size_t i = 0; i < result.domainSize(); ++i)
            result[i] = (*this)(i);
        return result;
    }

    // Get Frobenius/von Mises/max stress
    ScalarField<Real> sqrtStressMeasure() const {
        ScalarField<Real> result(microStress.domainSize());
        for (size_t i = 0; i < result.domainSize(); ++i)
            result[i] = sqrt((*this)(i));
        return result;
    }

    size_t size() const {
        return microStress.domainSize();;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Public data members
    ////////////////////////////////////////////////////////////////////////////
    MinorSymmetricRank4Tensor<N> Cbase;
    std::vector<Stress> stressField;
};

template<size_t N, bool _majorSymmCBase>
MicroscopicStress<N> MicroscopicFrobeniusStress(
        const ElasticityTensor<Real, N, _majorSymmCBase> &Cbase, // Allow non major-symmetric to handle the von Mises case
        std::vector<Stress> &stressField)
{
    MicroscopicStress<N> result;
    result.stressField = stressField;
    return result;
}

// Global worst-case objective of the form
// int j(s, x) dV
// Currently, j (and s) are considered piecewise constant per-element.
template<size_t N, class Integrand>
struct IntegratedWorstCaseObjective {
    using SMF = typename WorstCaseStress<N>::SMF;

    IntegratedWorstCaseObjective() { }
    template<class Mesh> IntegratedWorstCaseObjective(const Mesh &m, const MicroscopicStress<N>  &stress) { setPointwiseStress(m, stress); }
    template<class Mesh> IntegratedWorstCaseObjective(const Mesh &m,       MicroscopicStress<N>  &stress) { setPointwiseStress(m, std::move(stress)); }

    template<class Mesh> void setPointwiseStress(const Mesh &m, const MicroscopicStress<N>  &stress) {
        microStress = stress;
        integrand.init(m, microStress);
        m_updateEvalCache(m);
    }

    template<class Mesh> void setPointwiseStress(const Mesh &m, MicroscopicStress<N> &stress) {
        microStress = std::move(stress);
        integrand.init(m, microStress);
        m_updateEvalCache(m);
    }

    // Evaluate objective by integrating over m
    Real evaluate() const { return m_evalCache; }

    // NOTE: these are *NOT* scaled to account for possible reflected base cell
    // copies (computes just the pointwise value).
    ScalarField<Real> integrandValues() const {
        const size_t numElems = microStress.size();
        ScalarField<Real> result(numElems);
        for (size_t i = 0; i < numElems; ++i)
            result(i) = integrand.j(microStress(i), i);
        return result;
    }

    // Sensitivity of global objective integrand to cell problem fluctuation
    // strains (rank 2 tensor field tau in the writeup).
    // tau = j_prime 2 * sigma : C^base
    // NOTE: this is *NOT* scaled to account for possible reflected base cell
    // copies (computes just the pointwise value).
    void compute_tau(size_t kl, SMF &result) const {
        microStress.sensitivityToCellStrain(result);
        for (size_t e = 0; e < microStress.size(); ++e)
            result(e) *= integrand.j_prime(microStress(e), e);
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
    adjointDeltaJ(const NonPeriodicCellOps<Sim> &nonPeriodicCellOps) const {
        // Dilation and delta strain terms
        const auto mesh = nonPeriodicCellOps.mesh();
        size_t nv = mesh.numVertices();

        BENCHMARK_START_TIMER_SECTION("Adjoint Cell Problem");
        BENCHMARK_START_TIMER("Compute Tau");
        // Cache tau and adjoint solution
        SMF tau;
        compute_tau(tau);
        BENCHMARK_STOP_TIMER("Compute Tau");

        // rho is same as p in formulas
        VectorField<Real, N> rho;
        rho = NonPeriodicOps.solveAdjointCellProblem(NonPeriodicOps.sim().perElementStressFieldLoad(tau));

        // get value for displacements
        VectorField<Real, N> &u = nonPeriodicCellOps.displacement();

        BENCHMARK_STOP_TIMER_SECTION("Adjoint Cell Problem");

        BENCHMARK_START_TIMER_SECTION("Dilation and Delta strain");

        using OF = ScalarOneForm<N>;
        OF delta_j(nv);
        delta_j.clear();

#if USE_TBB
        tbb::combinable<OF> sum(delta_j);
#endif

        auto accumElementContrib = [&](size_t ei) {
            auto e = mesh.element(ei);
#if USE_TBB
            OF &result = sum.local();
#else
            OF &result = delta_j;
#endif
            // j contribution to dilation integrand

            // Part 1: int_w [j - strain(p):stress] grad(lambda_m) dx
            Real dilationIntegrand = integrand.j(microStress(e.index()), e.index());
            const auto C = e->E();
            using Strain = typename Sim::Strain;
            Strain strain_rho, strain_u;
            e->strain(e, rho, strain_rho);
            e->strain(e,   u, strain_u);

            using Stress = typename Sim::Stress;
            Stress stress_rho, stress_u;
            for (size_t i = 0; i < Stress::numNodalValues; ++i) {
                stress_rho[i] = C.doubleContract(strain_rho[i]);
                stress_u  [i] = C.doubleContract(strain_u  [i]);
            }

            // contribution to dilation integrand.
            // Part 1b: - int_w strain(p):stress grad(lambda_m) dx
            dilationIntegrand -= Quadrature<N, 2 * Strain::Deg>::integrate(
                    [&] (const EvalPt<N> &pt) {
                        return strain_rho(pt).doubleContract(stress_u(pt));
                    });

            // Finish Part 1: complete integral of dilatationIntegrand
            for (auto v : e.vertices()) {
                result(v.index()) += dilationIntegrand * e->volume()
                                     * e->gradBarycentric().col(v.localIndex());
            }

            // Delta strain term
            // Part 2: int_w [grad(lambda_m) . (stress * p_n + (strain(p):C - tau) u)] grad(phi_n) dx
            for (auto n : e.nodes()) {
                auto gradPhi_n = e->gradPhi(n.localIndex());

                // glam functional: tau * u - (stress * p_n + strain(p):C * u)
                Interpolant<VectorND<N>, N, Stress::Deg> glam_functional;
                glam_functional = tau(e.index()).contract(u(n.index()));
                for (size_t i = 0; i < glam_functional.size(); ++i) {
                    glam_functional[i] -= stress_rho[i].contract(  u(n.index()))
                                        + stress_u  [i].contract(rho(n.index()));
                }

                // compute grad of lambda_m and grad of phi_n and integrate:
                // -[grad(lambda_m) . glam_functional] grad(phi_n)
                for (auto v_m : e.vertices()) {
                    auto gradLam_m = e->gradBarycentric().col(v_m.localIndex());

                    // compute grad(lambda_m) . glam_functional
                    Interpolant<Real, N, Strain::Deg> mu_grad_lam_m;
                    for (size_t i = 0; i < mu_grad_lam_m.size(); ++i)
                        mu_grad_lam_m[i] = gradLam_m.dot(glam_functional[i]);

                    // finish computing Part 2 of integral and remove it from current result
                    result(v_m.index()) -= Quadrature<N, 2 * Strain::Deg>::integrate([&](const EvalPt<N> &pt) {
                        return (mu_grad_lam_m(pt) * gradPhi_n(pt)).eval();
                    }, e->volume());
                }
            }

        };

#if USE_TBB
        tbb::parallel_for(
                tbb::blocked_range<size_t>(0, mesh.numElements()),
                [&](const tbb::blocked_range<size_t> &r) {
                    for (size_t ei = r.begin(); ei < r.end(); ++ei) accumElementContrib(ei);
                });
        delta_j = sum.combine([](const OF &a, const OF &b) { return a + b; } );
#else
        for (auto e : sim.mesh().elements()) { accumElementContrib(e.index()); }
#endif

        BENCHMARK_STOP_TIMER_SECTION("Dilation and Delta strain");

        BENCHMARK_START_TIMER_SECTION("Surface traction term");

        // TODO: Is the surface integral term correct?
        // Part 3: int_(delta w) (rho . T) grad(lambda_m) dS
        for (auto be : m_mesh.boundaryElements()) {
            Real area = be->volume();
            Real boundaryIntegrand = pho.dot(be->neumannTraction);

            for (auto v : be.vertices()) {
                result(v.index()) += boundaryIntegrand * area * be->gradBarycentric().col(v.localIndex());
            }
        }

        BENCHMARK_STOP_TIMER_SECTION("Surface traction term");

        return delta_j;
    }

    ////////////////////////////////////////////////////////////////////////////
    // Public data members
    ////////////////////////////////////////////////////////////////////////////
    MicroscopicStress<N> microStress;
    Integrand integrand;

private:
    Real m_evalCache = 0.0;
    template<class Mesh> void m_updateEvalCache(const Mesh &m) {
        assert(m.numElements() == microStress.size());
        m_evalCache = 0;
        for (auto e : m.elements())
            m_evalCache += integrand.j(microStress(e.index()), e.index()) * e->volume();
    }
};

// int_omega micro_stress dV
// i.e. j(s, x) = s -> j' = 1.
struct MicroscopicStressIntegrandTotal {
    template<size_t N, class Mesh>
    void init(const Mesh &/* m */, const MicroscopicStress<N> &/* microStress */) { }

    // Derivative of global objective integrand wrt worst case stress.
    static Real j(Real microStress, size_t /* x_i */) { return microStress; }
    static Real j_prime(Real /* microStress */, size_t /* x_i */) { return 1; }
};

// int_omega worst_case_stress^p dV
// i.e. j(s, x) = s^p -> j' = p s^(p - 1).
struct MicroscopicStressIntegrandLp {
    template<size_t N, class Mesh>
    void init(const Mesh &/* m */, const MicroscopicCaseStress<N> &/* microStress */) { }

    MicroscopicStressIntegrandLp() { }

    Real j(Real microStress, size_t /* x_i */) const { return pow(microStress, p); }
    // Derivative of global objective integrand wrt stress.
    Real j_prime(Real microStress, size_t /* x_i */) const { return p * pow(microStress, p - 1); }
    Real p = 1.0;
};


// Take the pth root of a stress objective function.
// Based on chain rule, all derivatives are weighted by
// 1/p f^((1 / p) - 1), where f = SubObjective::eval();
template<class SubObjective>
struct PthRootObjective : public SubObjective {
    using Base = SubObjective;
    // Forward all constructor args to SubObjective
    template<typename... Args>
    PthRootObjective(Args&&... args)
            : SubObjective(std::forward<Args>(args)...) { }

    Real evaluate() const {
        if (p == 1) return SubObjective::evaluate();
        return pow(SubObjective::evaluate(), 1.0 / p);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Discrete shape derivative under mesh vertex perturbation delta_p
    // (forward mode).
    ////////////////////////////////////////////////////////////////////////////
    template<class Sim>
    Real deltaJ(const Sim &sim,
                const std::vector<typename Sim::VField> &w,
                const typename Sim::VField &delta_p) const {
        Real pder = SubObjective::deltaJ(sim, w, delta_p);
        if (p == 1) return pder;
        return m_gradientScale() * pder;
    }

    template<class Sim>
    ScalarOneForm<Sim::N>
    adjointDeltaJ(const NonPeriodicOperations<Sim> &NonPeriodicOps) const {
        ScalarOneForm<Sim::N> pder = SubObjective::adjointDeltaJ(NonPeriodicOps);
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

#endif //MICROSCOPICSTRESS_H
