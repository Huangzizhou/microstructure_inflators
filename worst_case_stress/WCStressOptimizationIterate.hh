////////////////////////////////////////////////////////////////////////////////
// WCStressOptimizationIterate.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Encapsulates the state of a worst-case stress pattern optimization
//      iterate and provides objective/gradient/etc.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/06/2015 23:16:57
////////////////////////////////////////////////////////////////////////////////
#ifndef WCSTRESSOPTIMIZATIONITERATE_HH
#define WCSTRESSOPTIMIZATIONITERATE_HH

#include <functional>
    
#include <PeriodicHomogenization.hh>
#include <MSHFieldWriter.hh>

#include "../pattern_optimization/PatternOptimizationIterate.hh"
#include "../pattern_optimization/ShapeVelocityInterpolator.hh"
#include "WorstCaseStress.hh"
#include "WCSObjective.hh"

#include "WCStressOptimizationConfig.hh"

#include <Laplacian.hh>
#include <MassMatrix.hh>

namespace WCStressOptimization {

// Implemented outside a member method so that a subclasses' static overrides
// of the evaluation/gradient computations will be called.
// (E.g. BoundaryPerturbationIterate must call its own gradp_JS since it
// overrides the default parameter shape velocity gradient approach).
template<class _WCSIterate>
void writeIterateDescription(std::ostream &os, _WCSIterate &it, bool printParams = false) {
    if (printParams) {
        auto params = it.getParams();
        os << "p:";
        for (size_t i = 0; i < params.size(); ++i)
            os << "\t" << params[i];
        os << std::endl;
    }

    os << "moduli:\t";
    it.elasticityTensor().printOrthotropic(os);
    os << "anisotropy:\t" << it.elasticityTensor().anisotropy() << std::endl;

    os << "printable:\t" << it.printable() << std::endl;

    BENCHMARK_START_TIMER("Evaluate JS/WCS");
    Real JS = it.evaluateJS(), WCS = it.evaluateWCS();
    os << "JS:\t"     << JS  << std::endl;
    os << "WCS:\t"    << WCS << std::endl;
    if (it.fullObjective().hasTargetVolume())
        os << "JVol:\t"   << it.evaluateJVol()  << std::endl;
    os << "Volume:\t" << it.mesh().volume() << std::endl;
    os << "Max Ptwise WCS:\t" << sqrt(it.wcsObjective().wcStress.stressMeasure().maxMag()) << std::endl;
    BENCHMARK_STOP_TIMER("Evaluate JS/WCS");

    auto   JS_p = it.gradp_JS();
    auto  WCS_p = it.gradientWCS_adjoint();
    ScalarField<Real> JVol_p(WCS_p.domainSize());
    JVol_p.clear();

    os << "||grad_p JS||:\t"  <<  JS_p.norm() << std::endl;
    os << "||grad_p WCS||:\t" << WCS_p.norm() << std::endl;

    if (it.fullObjective().hasTargetVolume()) {
        JVol_p = it.gradp_JVol();
        os << "||grad_p Jvol||:\t" << JVol_p.norm() << std::endl;
    }

    os << "Normalized WCS:\t" << WCS / it.fullObjective().initialWCS() << std::endl;
    os << "Normalized JS:\t" << it.evaluateJS() / it.fullObjective().targetSNormSq() << std::endl;
    if (it.fullObjective().hasTargetVolume())
        os << "Normalized JVol:\t" << it.evaluateJVol() / it.fullObjective().targetVolSq() << std::endl;

    // Output composite iterate stats.
    if (it.fullObjective().hasTargetVolume()) {
        os << "J_full:\t" << it.fullObjective().eval(JS, WCS, it.evaluateJVol()) << std::endl;
        os << "||grad_p J_full||:\t" << it.fullObjective().evalGradient(JS_p, WCS_p, JVol_p).norm() << std::endl;
    }
    else {
        os << "J_full:\t" << it.fullObjective().eval(JS, WCS) << std::endl;
        os << "||grad_p J_full||:\t" << it.fullObjective().evalGradient(JS_p, WCS_p).norm() << std::endl;
    }
    os << std::endl;
}

// _BypassParameterVelocity: hack to avoid computing the parameter velocity when
// there are too many parameters. This is only intended to be used by subclasses 
// that know how to compute objective gradients without explicitly forming the
// parameter velocity vectors (e.g. BoundaryPerturbationIterate)
template<class Sim, class WCSObjective = PthRootObjective<IntegratedWorstCaseObjective<Sim::N, WCStressIntegrandLp>>, bool _BypassParameterVelocity = false>
class Iterate : public PatternOptimization::Iterate<Sim, _BypassParameterVelocity> {
    using Base = PatternOptimization::Iterate<Sim, _BypassParameterVelocity>;
public:
    static constexpr size_t N = Sim::N;
    using ETensor = typename Sim::ETensor;
    using SField = ScalarField<Real>;

template<class _Inflator>
    // WARNING: updates fullObjective with init WCS on first construction (hack)
    Iterate(_Inflator &inflator, size_t nParams, const double *params,
            Objective<N> &fullObjective)
        : Base(inflator, nParams, params, fullObjective.targetS, true /* Always keep w_ij */ ),
          m_fullObjective(fullObjective)
    {
        if (m_fullObjective.hasTargetVolume())
            this->setTargetVolume(m_fullObjective.targetVolume());

        // Worst case stress currently assumes that the base material is
        // constant, so we can read it off a single element.
        m_wcs_objective.setPointwiseWCS(m_sim->mesh(),
            worstCaseFrobeniusStress(m_sim->mesh().element(0)->E(), Base::S,
                PeriodicHomogenization::macroStrainToMicroStrainTensors(w_ij, *m_sim)));

        if (!m_fullObjective.hasInitialWCS()) {
            m_fullObjective.setInitialWCS(evaluateWCS());
            std::cout << "Initial WCS: " << m_fullObjective.initialWCS() << std::endl;
        }

        // Get the shape velocity vector fields needed for the discrete shape derivative.
        if (!_BypassParameterVelocity)
            m_bdry_svels = inflator.shapeVelocities(m_sim->mesh());
    }

    // Evaluate the global worst case stress objective
    Real evaluateWCS() const { return m_wcs_objective.evaluate(); }

    Real evaluateJFull() const {
        if (m_fullObjective.hasTargetVolume()) return m_fullObjective.eval(Base::evaluateJS(), evaluateWCS(), Base::evaluateJVol());
        else                                   return m_fullObjective.eval(Base::evaluateJS(), evaluateWCS());
    }

    // WARNING: paste this into any subclass that overrides gradient evaluation.
    SField gradp_JFull() const {
        if (m_fullObjective.hasTargetVolume()) return m_fullObjective.evalGradient(Base::gradp_JS(), gradientWCS_adjoint(), Base::gradp_JVol());
        else                                   return m_fullObjective.evalGradient(Base::gradp_JS(), gradientWCS_adjoint());
    }

    // Direct differentation version of gradient
    ScalarField<Real> gradientWCS_direct() const {
        ScalarField<Real> result(Base::m_params.size());
        for (size_t p = 0; p < result.domainSize(); ++p)
            result[p] = m_wcs_objective.directDerivative(*m_sim, w_ij, m_vn_p[p]);
        return result;
    }

    // Direct differentation version of gradient
    Real gradientWCS_direct_component(size_t p) const {
        return m_wcs_objective.directDerivative(*m_sim, w_ij, m_vn_p.at(p));
    }

    // Adjoint method version of gradient: evaluates linear functional form of
    // shape derivative
    ScalarField<Real> gradientWCS_adjoint() const {
        ScalarField<Real> result(Base::m_params.size());

        auto shapeDerivative = m_wcs_objective.gradient(*m_sim, w_ij);

        using  SDInterp = typename std::decay<decltype(shapeDerivative[0])>::type;
        using NSVInterp = typename std::decay<decltype(      m_vn_p[0][0])>::type;
        static_assert(SDInterp::K == NSVInterp::K,
                "Invalid boundary interpolant simplex dimension");

        // Boundary interpolant...
        for (size_t p = 0; p < result.domainSize(); ++p) {
            result[p] = 0;
            for (auto be : m_sim->mesh().boundaryElements()) {
                const auto &vn = m_vn_p[p][be.index()];
                const auto &sd = shapeDerivative[be.index()];
                result[p] += Quadrature<SDInterp::K, SDInterp::Deg + NSVInterp::Deg>::
                    integrate([&] (const VectorND<be.numVertices()> &pt) {
                        return vn(pt) * sd(pt);
                    }, be->volume());
            }
        }

        return result;
    }

    Real gradientWCS_discrete_forward(size_t p) const {
        assert(!_BypassParameterVelocity);
        // Interpolate boundary shape velocity into the interior vertices to
        // improve gradient accuracy. (TODO: Make this an option).

        // MSHBoundaryFieldWriter bdryWriter("bdry_svel_debug.msh", m_sim->mesh());
        // bdryWriter.addField("bdry_svel", m_bdry_svels.at(p));

        ShapeVelocityInterpolator interpolator(*m_sim);
        auto svel = interpolator.interpolate(*m_sim, m_bdry_svels.at(p));

        // MSHFieldWriter writer("svel_debug.msh", m_sim->mesh());
        // writer.addField("svel", svel);

        return m_wcs_objective.deltaJ(*m_sim, w_ij, svel);
    }

    Real gradientWCS_discrete_adjoint(size_t p) const {
        assert(!_BypassParameterVelocity);

        auto delta_j    = m_wcs_objective.adjointDeltaJ(*m_sim, w_ij);
        ShapeVelocityInterpolator interpolator(*m_sim);
        auto delta_j_vb = interpolator.adjoint(*m_sim, delta_j);

        // MSHBoundaryFieldWriter bdryWriter("bdry_adjoint_debug.msh", m_sim->mesh());
        // bdryWriter.addField("delta_j_vb", delta_j_vb, DomainType::PER_NODE);

        const auto &bsvel = m_bdry_svels.at(p);
        size_t numBV = m_sim->mesh().numBoundaryVertices();
        assert(delta_j_vb.domainSize() == numBV);
        assert(     bsvel.domainSize() == numBV);
        Real dJ = 0;
        for (size_t bvi = 0; bvi < numBV; ++bvi)
            dJ += delta_j_vb(bvi).dot(bsvel(bvi));

#if 0
        // Debug steepest descent directions g_vol and g_bdry,
        // by ensuring <g_vol, v_vol> = <g_bdry, v_bdry> = dJ_bdry[v_bdry].
        Real vol_representative_dJ = 0;
        auto svel = interpolator.interpolate(*m_sim, m_bdry_svels.at(p));
        auto grad = steepestDescentVolumeVelocity();
        for (auto e : m_sim->mesh().elements()) {
            Interpolant<VectorND<N>, N, 1> v_p, g;
            for (auto v : e.vertices()) {
                v_p[v.localIndex()] = svel(v.index());
                g  [v.localIndex()] = grad(v.index());
            }

            vol_representative_dJ += Quadrature<N, 2>::integrate([&](const VectorND<e.numVertices()> &pt) {
                    return v_p(pt).dot(g(pt));
                }, e->volume());
        }

        Real bdry_representative_dJ = 0;
        auto bdry_grad = steepestDescentBoundaryVelocity();
        for (auto be : m_sim->mesh().boundaryElements()) {
            Interpolant<VectorND<N>, N - 1, 1> v_p, g;
            for (auto bv : be.vertices()) {
                v_p[bv.localIndex()] = bsvel(bv.index());
                g  [bv.localIndex()] = bdry_grad(bv.index());
            }
            bdry_representative_dJ += Quadrature<N - 1, 2>::integrate([&](const VectorND<be.numVertices()> &pt) {
                    return v_p(pt).dot(g(pt));
                }, be->volume());
        }

        std::cout << std::endl;
        std::cout << "correct dJ: " << dJ << std::endl;
        std::cout << "Riesz representative dJ (vol): "  <<  vol_representative_dJ << std::endl;
        std::cout << "Riesz representative dJ (bdry): " << bdry_representative_dJ << std::endl;
        std::cout << std::endl;
#endif

        return dJ;
    }

    // Compute from the objective's differential one-form a mesh-independent
    // steepest-descent volume velocity field ("steepest" with respect to the
    // period cell geometry's L^2 norm).
    // This amounts to computing a periodic velocity field, g, whose integrated
    // inner product over the mesh against velocity v gives the change in
    // objective (i.e. a "Riesz representative" for dJ):
    //      int_omega g(x) . v(x) dx := dJ[v] = - sum_i delta_j[i] . v[i]
    // where v(x) := v[k] * phi^k(x),   (periodic)
    //       g(x) := g[k] * phi^k(x),   (periodic)
    // and phi_k is a scalar vertex shape function.
    //
    // Correctly handling the periodicity requires care:
    // Notice the sum over i in the one-form evaluation is over **all**
    // vertices, so it's effectively summing all partial contributions in
    // delta_j[i] for periodically identified vertices.
    //
    // Introducing "selection" matrix S that distributes the reduced periodic DoF
    // values to the corresponding vertices (v = S v_dof):
    //   dJ[v] = sum_i delta_j[i] . (S v_dof)[i] = sum_i (S^T delta_j)[i] . v_dof[i]
    // (here, S^T sums the partial contributions each identified vertex into a
    //  per-dof total).
    // Likewise, we can write v(x) = (S v_dof)[k] phi^k(x),
    //                        g(x) = (S g_dof)[k] phi^k(x) ==>
    //      int_omega g(x) . v(x) dx := -dJ[v] = -sum_i delta_j[i] . v[i]
    //    = g_dof[k] (int_omega S_sk phi^s(x) phi^t(x) S_tl) . v_dof[l]
    //   := g_dof[k] . v_dof[l] [M_dof]_kl
    // with reduced mass matrix M_dof = S^T M S.
    // Finally, for this inner product to equal -dJ[v] for all v, we see:
    //      sum_i (M_dof g_dof)[i] . v_dof[i] = -sum_i (S^T delta_j)[i] . v_dof[i]
    //  ==> g_dof = -M_dof^{-1} S^T delta_j
    // and we can recover the velocity field coefficients g = S g_dof
    typename Sim::VField steepestDescentVolumeVelocity() const {
        auto delta_j = m_wcs_objective.adjointDeltaJ(*m_sim, w_ij);
        // Determine S, the map from reduced periodic DoFs to vertices
        // Unfortunately, sim's periodic boundary conditions are on nodes, not
        // just vertices, so we have to reindex to get contiguous variables
        const auto &mesh = m_sim->mesh();
        size_t nv = mesh.numVertices();
        constexpr size_t NO_VAR = std::numeric_limits<size_t>::max();
        std::vector<size_t> dofForSimNodalDoF(m_sim->numDoFs(), NO_VAR);
        std::vector<size_t> dofForVertex;
        dofForVertex.reserve(nv);
        size_t numDoFs = 0;
        for (auto v : mesh.vertices()) {
            size_t &dof = dofForSimNodalDoF.at(m_sim->DoF(v.node().index()));
            if (dof == NO_VAR) dof = numDoFs++;
            dofForVertex.push_back(dof);
        }

        auto M = MassMatrix::construct<1>(mesh);
        // ==> S M S^T
        M.reindexVariables(numDoFs, dofForVertex);
        SPSDSystem<Real> M_dof(M);

        typename Sim::VField g(nv);
        for (size_t c = 0; c < N; ++c) {
            // Apply -S^T
            std::vector<Real> neg_S_t_delta_j(numDoFs, 0.0);
            for (auto v : mesh.vertices())
                neg_S_t_delta_j.at(dofForVertex[v.index()]) -= delta_j(v.index())[c];

            auto g_dof = M_dof.solve(neg_S_t_delta_j);

            // Apply S 
            for (auto v : mesh.vertices())
                g(v.index())[c] = g_dof.at(dofForVertex.at(v.index()));
        }
        
        return g;
    }

    // Compute from the objective's differential one-form a mesh-independent
    // steepest-descent boundary velocity field ("steepest" with respect to the
    // boundary geometry's L^2 norm):
    //      min dJ[g]                   (dJ[g] = sum_i delta_j_bdry[i] . g[i])
    //       g         ==> delta_j_bdry + 2 l M g = 0
    //   g^T M g = h^2    
    //      ==> g \propto -M^-1 delta_j_bdry
    // Here M is the boundary mass matrix (g^T M g = int_bdry g(x) . g(x) dA(x),
    // where g(x) is the interpolated scalar field corresponding to vector g).
    // This g can also be interpreted as the "Riesz representative" for dJ.
    // Periodicity details:
    //    We seek a periodic g, and the boundary integral defining M is over
    //    the "true" boundary (i.e. the boundary remaining after stitching
    //    together identified vertices). We simplify the problem by reducing to
    //    "periodic dof" variables.
    //    Introducing matrix S that selects the dof value for each vertex:
    //          g = S g_dof
    //    (Note that S^T sums values over all identified vertices and places
    //     them on the corresponding periodic dof.)
    //
    //    delta_j_bdry is non-periodic; it contains partial contributions on
    //    each vertex (so it correctly evaluates the one-form when dotted with
    //    a *periodic* discrete boundary nodal vector field). Thus the summed
    //    quantity S^T delta_j_bdry gives the one-form acting on periodic
    //    boundary node displacements.
    //    
    //    Letting M_dof = S^T M S, we arrive at the equation:
    //         g_dof \propto -  M_dof^-1 S^T delta_j_bdry
    //    ==>  g     \propto -S M_dof^-1 S^T delta_j_bdry
    // We let the caller decide on a normalization by simply computing:
    // @return      -S M_dof^-1 S^T delta_j_bdry
    typename Sim::VField steepestDescentBoundaryVelocity() const {
        ShapeVelocityInterpolator interpolator(*m_sim);
        auto delta_j      = m_wcs_objective.adjointDeltaJ(*m_sim, w_ij);
        auto delta_j_bdry = interpolator.adjoint(*m_sim, delta_j);
        
        // Extract boundary vertex periodicity from m_sim.
        // Unfortunately, sim's periodic boundary conditions are on volume
        // nodes, not boundary vertices, so we have to reindex to get
        // contiguous variables
        const auto mesh = m_sim->mesh();
        size_t nbv = mesh.numBoundaryVertices();
        std::vector<size_t> dofForBdryVertex; dofForBdryVertex.reserve(nbv);

        constexpr size_t NO_VAR = std::numeric_limits<size_t>::max();
        std::vector<size_t> dofForSimNodalDoF(m_sim->numDoFs(), NO_VAR);
        size_t numDoFs = 0;
        for (auto bv : mesh.boundaryVertices()) {
            size_t simDoF = m_sim->DoF(bv.volumeVertex().node().index());
            size_t &dof = dofForSimNodalDoF.at(simDoF);
            if (dof == NO_VAR)
                dof = numDoFs++;
            dofForBdryVertex.push_back(dof);
        }

        // Determine which variables correspond to true boundary vertices
        // (if any non-periodic boundary element touches them), and extract an
        // array for testing periodic boundary element membership.
        std::vector<bool> isTrueBoundaryVertex(mesh.numBoundaryVertices(), false),
                          isPeriodicBE(mesh.numBoundaryElements(), false);
        for (auto be : mesh.boundaryElements()) {
            isPeriodicBE[be.index()] = be->isPeriodic;
            if (be->isPeriodic) continue;
            for (size_t i = 0; i < be.numVertices(); ++i)
                isTrueBoundaryVertex.at(be.vertex(i).index()) = true;
        }

        // Determine dofs corresponding to false boundary vertices (to be constrained)
        std::vector<bool> dofVisited(numDoFs, false); // Add each dof only once
        std::vector<size_t> falseBoundaryDoFs;
        for (auto bv : mesh.boundaryVertices()) {
            if (isTrueBoundaryVertex.at(bv.index())) continue;

            size_t dof = dofForBdryVertex.at(bv.index());
            if (dofVisited.at(dof)) continue;
            falseBoundaryDoFs.push_back(dof);
            dofVisited[dof] = true;
        }

        // Build M_dof, ignoring contributions from periodic boundary elements
        auto M = MassMatrix::construct<1>(mesh.boundary(), false, isPeriodicBE);
        // ==> S M S^T
        M.reindexVariables(numDoFs, dofForBdryVertex);
        SPSDSystem<Real> M_dof(M);

        // Enforce zero velocity on false boundary vertices
        // (Since we skipped periodic bdry elements' contributions, the system
        //  would be singular otherwise).
        M_dof.fixVariables(falseBoundaryDoFs,
                           std::vector<Real>(falseBoundaryDoFs.size(), 0.0));

        typename Sim::VField g(nbv);
        for (size_t c = 0; c < N; ++c) {
            // Apply -S^T
            std::vector<Real> neg_S_t_delta_j_bdry(numDoFs, 0.0);
            for (auto bv : mesh.boundaryVertices()) {
                neg_S_t_delta_j_bdry.at(dofForBdryVertex.at(bv.index()))
                    -= delta_j_bdry(bv.index())[c];
            }
            
            auto g_dof = M_dof.solve(neg_S_t_delta_j_bdry);

            // Apply S
            for (auto bv : mesh.boundaryVertices())
                g(bv.index())[c] = g_dof.at(dofForBdryVertex.at(bv.index()));
        }

        return g;
    }

    // WARNING: paste this into any subclass that overrides gradient/objective evaluation.
    void writeDescription(std::ostream &os) const {
        writeIterateDescription(os, *this, true);
    }

    // Also output intermediates of derivative computation wrt parameter
    // "derivativeComponent" if requested.
    void writeMeshAndFields(const std::string &name, int derivativeComponent = -1, bool fullDegreeFieldOutput = false) const {
        MSHFieldWriter writer(name, m_sim->mesh(), /* linearSubsample: */ !fullDegreeFieldOutput);
        // typename Sim::VField outField;
        // for (size_t kl = 0; kl < flatLen(N); ++kl) {
        //     // Subtract off average displacements so that fields are comparable
        //     // across meshes.
        //     outField = w_ij[kl];
        //     outField -= outField.mean();
        //     writer.addField("w " + std::to_string(kl), outField);
        //     writer.addField("we " + std::to_string(kl), m_sim->averageStrainField(w_ij[kl]));
        // }
        ScalarField<Real> j = m_wcs_objective.integrandValues();

        writer.addField("Pointwise WCS", m_wcs_objective.wcStress.sqrtStressMeasure());
        writer.addField("j", j);

        writer.addField("Steepest Descent VVel", steepestDescentVolumeVelocity(), DomainType::PER_NODE);
        auto bdryVel = steepestDescentBoundaryVelocity();
        VectorField<Real, N> xferBdryVel(m_sim->mesh().numVertices());
        xferBdryVel.clear();
        for (auto v : m_sim->mesh().vertices()) {
            auto bv = v.boundaryVertex();
            if (!bv) continue;
            xferBdryVel(v.index()) = bdryVel(bv.index());
        }
        writer.addField("Steepest Descent BVel", xferBdryVel, DomainType::PER_NODE);

        if (derivativeComponent >= 0) {
#if 0
            std::vector<VectorField<Real, N>> dot_w;
            PeriodicHomogenization::fluctuationDisplacementShapeDerivatives(*m_sim, w_ij, m_vn_p[derivativeComponent], dot_w,
                WCStressOptimization::Config::get().projectOutNormalStress);
            typename WCSObjective::SMF tau;

            // Output full-degree tensor fields. (Wasteful since
            // strain fields are of degree - 1, but Gmsh/MSHFieldWriter
            // only supports full-degree ElementNodeData).
            using UpsampledTensorInterp = SymmetricMatrixInterpolant<typename Sim::SMatrix, Sim::N, Sim::Degree>;
            using OrigTensorInterp = typename Sim::Strain;
            auto upsampledField = [&](const std::vector<OrigTensorInterp> &orig) -> std::vector<UpsampledTensorInterp> {
                std::vector<UpsampledTensorInterp> upsampledField; upsampledField.reserve(orig.size());
                for (const auto s: orig) upsampledField.emplace_back(s);
                return upsampledField;
            };

            for (size_t i = 0; i < dot_w.size(); ++i) {
                m_wcs_objective.tau_kl(i, tau);

                // writer.addField("w_ij " + std::to_string(i), w_ij[i]);
                // writer.addField("dw_ij " + std::to_string(i), dot_w[i]);

                if (fullDegreeFieldOutput && (Sim::Degree > 1)) {
                    writer.addField("strain w_ij "  + std::to_string(i), upsampledField(m_sim->strainField(w_ij[i])));
                    writer.addField("strain dw_ij " + std::to_string(i), upsampledField(m_sim->strainField(dot_w[i])));
                }
                else {
                    writer.addField("strain w_ij "  + std::to_string(i), m_sim->averageStrainField(w_ij[i]));
                    writer.addField("strain dw_ij " + std::to_string(i), m_sim->averageStrainField(dot_w[i]));
                }
                writer.addField("tau_kl " + std::to_string(i), tau);
            }
            SField dj = m_wcs_objective.directIntegrandDerivative(*m_sim, w_ij, dot_w, m_vn_p[derivativeComponent]);
            writer.addField("dj", dj);
#endif
        }
    }
    
    const WCSObjective &wcsObjective() const { return m_wcs_objective; }
    const WCStressOptimization::Objective<N> &fullObjective() const { return m_fullObjective; }

protected:
    WCStressOptimization::Objective<N> &m_fullObjective;
    WCSObjective m_wcs_objective;
    // Full, vector-valued shape velocity boundary vector fields (if param
    // velocity isn't bypassed)
    std::vector<VectorField<Real, N>> m_bdry_svels;
    using Base::m_sim;
    using Base::w_ij;
    using Base::m_vn_p;
};

}

#endif /* end of include guard: WCSTRESSOPTIMIZATIONITERATE_HH */
