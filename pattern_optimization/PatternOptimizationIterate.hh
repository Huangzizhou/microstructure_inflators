////////////////////////////////////////////////////////////////////////////////
// PatternOptimizationIterate.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Encapsulates the state of a pattern optimization iterate and provides
//      objective/gradient/etc.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/26/2014 19:04:20
////////////////////////////////////////////////////////////////////////////////
#ifndef PATTERNOPTIMIZATIONITERATE_HH
#define PATTERNOPTIMIZATIONITERATE_HH

#include <EdgeFields.hh>
#include <MSHFieldWriter.hh>

#include <iostream>
#include <cstdio>
#include <cassert>
#include <memory>
#include <iostream>
#include <iomanip>
#include <memory>
#include <tuple>

#include <MeshIO.hh>
#include <TriMesh.hh>
#include <filters/subdivide.hh>

#include "PatternOptimizationConfig.hh"

#include <boost/optional.hpp>
#include <Future.hh>

namespace PatternOptimization {

////////////////////////////////////////////////////////////////////////////////
// Geometry post-processing stage. Allow modification of the inflation geometry
// before building the simulator (but without changing the inflator's copy of
// the geometry. This hack enables support, for e.g., 1->4 subdivision in the
// 2D case with zero overhead when not requested.
////////////////////////////////////////////////////////////////////////////////
template<size_t N>
struct GeometryPostProcessor {
    template<class _Inflator>
    static
    std::tuple<std::unique_ptr<std::vector<MeshIO::IOVertex >>,
               std::unique_ptr<std::vector<MeshIO::IOElement>>> run(const _Inflator &/* i */) {
        return std::make_tuple(std::unique_ptr<std::vector<MeshIO::IOVertex >>(),
                               std::unique_ptr<std::vector<MeshIO::IOElement>>());
    }
};

template<>
struct GeometryPostProcessor<2> {
    template<class _Inflator>
    static
    std::tuple<std::unique_ptr<std::vector<MeshIO::IOVertex >>,
               std::unique_ptr<std::vector<MeshIO::IOElement>>> run(const _Inflator &i) {
        const auto &config = Config::get();
        if (config.fem2DSubdivRounds == 0) {
            return std::make_tuple(std::unique_ptr<std::vector<MeshIO::IOVertex >>(),
                                   std::unique_ptr<std::vector<MeshIO::IOElement>>());
        }
        else {
            // NOTE: subdividing here will break the shape derivative/gradient
            // comptuation (iterate will have a higher resolution boundary
            // discretization than the inflator, where parameter normal velocities are
            // defined).
            // TODO: implement subdivision for the boundary perturbation
            // inflator in a "reduced boundary perturbation" class instead.
            auto subv = Future::make_unique<std::vector<MeshIO::IOVertex >>(i.vertices());
            auto sube = Future::make_unique<std::vector<MeshIO::IOElement>>(i.elements());
            for (size_t i = 0; i < config.fem2DSubdivRounds; ++i) {
                ::TriMesh<SubdivVertexData<2>, SubdivHalfedgeData> m(*sube, subv->size());
                for (size_t vi = 0; vi < subv->size(); ++vi)
                    m.vertex(vi)->p = (*subv)[vi];
                subdivide(m, *subv, *sube);
            }
            return std::make_tuple(std::move(subv), std::move(sube));
        }
    }
};

// _BypassParameterVelocity: hack to avoid computing the parameter velocity when
// there are too many parameters. This is only intended to be used by subclasses 
// that know how to compute objective gradients without explicitly forming the
// parameter velocity vectors (e.g. BoundaryPerturbationIterate)
template<class _Sim, bool _BypassParameterVelocity = false>
struct Iterate {
    typedef typename _Sim::VField VField;
    typedef ScalarField<Real> SField;
    typedef typename _Sim::ETensor _ETensor;
    static constexpr size_t _N = _Sim::N;
    // 2 * (deg - 1) boundary element interpolant of an elasticity tensor
    // field used to represent the per-boundary-element value of the shape
    // derivative of the homogenized tensors.
    typedef PeriodicHomogenization::
            BEHTensorGradInterpolant<_Sim>          BEGradTensorInterpolant;
    // 2 * (deg - 1) boundary element interpolant of a scalar field used to
    // represent the per-boundary-element value of the shape derivative of a
    // scalar function.
    typedef Interpolant<Real, BEGradTensorInterpolant::K,
                    BEGradTensorInterpolant::Deg>   BEGradInterpolant;

    template<class _Inflator>
    Iterate(_Inflator &inflator, size_t nParams, const double *params,
            const _ETensor &targetS, bool keepFluctuationDisplacements = false)
        : m_targetS(targetS)
    {
        m_params.resize(nParams);
        for (size_t i = 0; i < nParams; ++i)
            m_params[i] = params[i];
        m_printable = inflator.isPrintable(m_params);

        // std::cout << "Inflating" << std::endl;
        BENCHMARK_START_TIMER_SECTION("Inflate");
        try {
            inflator.inflate(m_params);
        }
        catch (...) {
            // Hack to correct timer behavior--should probably use RAII
            BENCHMARK_STOP_TIMER_SECTION("Inflate");
            throw;
        }
        BENCHMARK_STOP_TIMER_SECTION("Inflate");
        // std::cout << "Inflated" << std::endl;

        // std::cout << "Checking geometry" << std::endl;
        if ((inflator.elements().size() == 0) || (inflator.vertices().size() == 0)) {
            std::cerr << std::setprecision(20);
            std::cerr << "Exception while inflating parameters" << std::endl;
            for (size_t i = 0; i < m_params.size(); ++i) std::cerr << m_params[i] << "\t";
            std::cerr << std::endl;
            throw std::runtime_error("Empty inflated geometry. Elements: "
                    + std::to_string(inflator.elements().size()) + ", Vertices: "
                    + std::to_string(inflator.vertices().size()));
        }
        // std::cout << std::endl;


        // std::cout << "Building Simulator" << std::endl;
        BENCHMARK_START_TIMER_SECTION("Eval");

        {
            std::unique_ptr<std::vector<MeshIO::IOVertex >> verts;
            std::unique_ptr<std::vector<MeshIO::IOElement>> elems;
            std::tie(verts, elems) = GeometryPostProcessor<_Sim::K>::run(inflator);
            if (elems && verts)
                m_sim = Future::make_unique<_Sim>(*elems, *verts);
            else {
                assert(!elems && !verts);
                m_sim = Future::make_unique<_Sim>(inflator.elements(),
                                                  inflator.vertices());
            }
        }
        // std::cout << "Done" << std::endl;
        // std::cout << "Homogenizing" << std::endl;

        try {
            PeriodicHomogenization::solveCellProblems(w_ij, *m_sim);
        }
        catch(std::exception &e) {
            std::cerr << "Cell problem solve failed: " << e.what() << std::endl;
            MeshIO::save("debug.msh", mesh());
            std::cerr << "Wrote geometry to 'debug.msh'" << std::endl;
            std::cerr << std::setprecision(19) << std::endl;
            std::cerr << "params:";
            for (size_t i = 0; i < m_params.size(); ++i) std::cerr << "\t" << m_params[i];
            std::cerr << std::endl;
            exit(-1);
        }

        // Shape velocities must be computed after periodic boundary conditions
        // are applied (i.e. solveCellProblems call)!
        if (!_BypassParameterVelocity)
            m_vn_p = inflator.computeShapeNormalVelocities(mesh());

        C = PeriodicHomogenization::homogenizedElasticityTensorDisplacementForm(w_ij, *m_sim);
        S = C.inverse();
        std::vector<BEGradTensorInterpolant> gradEh =
            PeriodicHomogenization::homogenizedElasticityTensorGradient(w_ij, *m_sim);
        // Compute compliance tensor gradient from elasticity tensor
        // gradient via chain rule
        // (doc/pattern_optimization/shape_derivative)
        m_gradS.resize(gradEh.size());
        for (size_t i = 0; i < gradEh.size(); ++i) {
            const auto &GE = gradEh[i];
                  auto &GS = m_gradS[i];
            // Compute each nodal value of the interpolant.
            for (size_t n = 0; n < GE.size(); ++n)
                GS[n] = -S.doubleDoubleContract(GE[n]);
        }

        if (!_BypassParameterVelocity) {
            // Precompute gradient of the compliance tensor
            m_gradp_S.resize(nParams); // Fill with zero tensors.
            for (size_t p = 0; p < nParams; ++p) {
                for (size_t bei = 0; bei < mesh().numBoundaryElements(); ++bei) {
                    auto be = mesh().boundaryElement(bei);
                    const auto &vn = m_vn_p[p][bei]; // parameter normal shape velocity interpolant
                    const auto &grad = m_gradS[bei];
                    m_gradp_S[p] += Quadrature<_Sim::K - 1, 1 + BEGradTensorInterpolant::Deg>::
                        integrate([&] (const VectorND<be.numVertices()> &pt) {
                            return vn(pt) * grad(pt);
                        }, be->volume());
                }
            }
        }

        m_diffS = S - m_targetS;

        // std::cout << "Done" << std::endl;

        if (!keepFluctuationDisplacements) w_ij.clear();

        BENCHMARK_STOP_TIMER_SECTION("Eval");
    }

    // To be used only if inflation fails--use a linear extrapolation to
    // estimate the objective/gradient at params
    void estimatePoint(size_t nParams, const double *params) {
        assert(nParams == m_params.size());
        m_estimateObjectiveWithDeltaP.resize(nParams);

        std::cerr << "WARNING, USING APPROXIMATE OBJECTIVE/GRADIENT AT DIST:";
        for (size_t p = 0; p < nParams; ++p) {
            m_estimateObjectiveWithDeltaP[p] = params[p] - m_params[p];
            std::cerr << "\t" << m_estimateObjectiveWithDeltaP[p];
        }
        std::cerr << std::endl;
    }
    void disableEstimation() { m_estimateObjectiveWithDeltaP.clear(); }

    // Evaluate compliance frobenius norm objective.
    Real evaluateJS() const {
        // Note: following is equivalent when computing the exact objective
        // (i.e. when m_estimateObjectiveWithDeltaP == {}):
        //      return 0.5 * m_diffS.quadrupleContract(m_diffS);
        // But linearly approximating the least squares objective is different
        // from linearly approximating the residual, and we prefer to do the
        // latter. 
        
        // Residual version
        Real result = 0;
        for (size_t i = 0; i < flatLen(_N); ++i) {
            for (size_t j = i; j < flatLen(_N); ++j) {
                Real r = residual(i, j);
                result += r * r;
            }
        }
        return 0.5 * result;
    }

    Real evaluateRT(SField &initialParams, double regularizationWeight) const {
    	Real result = 0.0;
		for (size_t p = 0; p < m_params.size(); ++p)
    		result += regularizationWeight * (m_params[p] - initialParams[p]) * (m_params[p] - initialParams[p]);
    	return result;
	}

    // S_ijkl - target_ijkl
    // EXCEPT when ignoreShear = true, in which case the "shear modulus
    // components" are zeroed out. (rows/cols >= _N)
    _ETensor diffS() const {
        if (PatternOptimization::Config::get().ignoreShear) {
            _ETensor zeroedShear = m_diffS;
            for (size_t i = _N; i < flatLen(_N); ++i)
                for (size_t j = i; j < flatLen(_N); ++j)
                    zeroedShear.D(i, j) = 0.0;
            return zeroedShear;
        }
        return m_diffS;
    }

    ////////////////////////////////////////////////////////////////////////
    /*! Computes grad(1/2 sum_ijkl (S_ijkl - target_ijlk|)^2) =
    //      (S_ijkl - target_ijlk) * grad(S_ikjl))
    //  @param[in]  target  S^* (target compliance tensor)
    //  @return Per-boundary-edge piecewise 2 * (deg - 1) scalar field giving
    //          steepest ascent normal velocity perturbation for JS
    *///////////////////////////////////////////////////////////////////////
    std::vector<BEGradInterpolant> shapeDerivativeJS() const {
        std::vector<BEGradInterpolant> grad(m_gradS.size());
        auto deltaS = diffS();
        for (size_t be = 0; be < m_gradS.size(); ++be) {
            // Compute each nodal value of the interpolant.
            const auto &GS = m_gradS[be];
                  auto &g = grad[be];
            for (size_t n = 0; n < GS.size(); ++n)
                g[n] = deltaS.quadrupleContract(GS[n]);
        }

        return grad;
    }

    // Computes grad_p(1/2 sum_ijkl (S_ijkl - target_ijlk|)^2) =
    //      (S_ijkl - target_ijlk) * grad_p(S_ikjl))
    SField gradp_JS() const {
        auto deltaS = diffS();
        SField result(m_params.size());
        for (size_t p = 0; p < m_params.size(); ++p)
            result[p] = deltaS.quadrupleContract(m_gradp_S[p]);
        return result;
    }

    SField gradp_RT(SField &initialParams, double regularizationWeight) const {
    	SField result(m_params.size());
    	for (size_t p = 0; p < m_params.size(); ++p)
    		result[p] = 2.0 * regularizationWeight * (m_params[p] - initialParams[p]);
    	return result;
	}

	SField getParams() const {
		SField result(m_params.size());
		for (size_t p = 0; p < m_params.size(); ++p)
			result[p] = m_params[p];
		return result;
	}


	// The (ij, kl)th residual (kl >= ij) for the nonlinear least squares (a
    // single term of the Frobenius distance). The terms are weighted so
    // that the squared norm of the residual vector corresponds to the
    // Frobenius norm of the rank 4 tensor difference S - S^*.
    Real residual(size_t ij, size_t kl) const {
        assert(kl >= ij);
        Real weight = 1.0;
        if (kl != ij) weight *= sqrt(2); // Account for lower triangle
        if (PatternOptimization::Config::get().ignoreShear) {
            if (ij >= _N) weight = 0.0; // Zero out shear components
            if (kl >= _N) weight = 0.0; // Zero out shear components
        }
        else {
            if (ij >= _N) weight *= sqrt(2); // Left shear doubler
            if (kl >= _N) weight *= sqrt(2); // Right shear doubler
        }
        Real result = weight * m_diffS.D(ij, kl);

        if (m_estimateObjectiveWithDeltaP.size() == m_params.size()) {
            for (size_t p = 0; p < m_params.size(); ++p)
                result += jacobian(ij, kl, p) * m_estimateObjectiveWithDeltaP[p];
        }
        return result;
    }

    // Derivative of residual(ij, kl) wrt parameter p:
    // d/dp (S_ijkl - target_ijkl) = d/dp S_ijkl = <gradS_ijkl, vn_p>
    // The terms are weighted in accordance with the residual weighting above.
    Real jacobian(size_t ij, size_t kl, size_t p) const {
        assert(kl >= ij);
        Real weight = 1.0;
        if (kl != ij) weight *= sqrt(2); // Account for lower triangle
        if (PatternOptimization::Config::get().ignoreShear) {
            if (ij >= _N) weight = 0.0; // Zero out shear components
            if (kl >= _N) weight = 0.0; // Zero out shear components
        }
        else {
            if (ij >= _N) weight *= sqrt(2); // Left shear doubler
            if (kl >= _N) weight *= sqrt(2); // Right shear doubler
        }
        return weight * m_gradp_S[p].D(ij, kl);
    }

    // Boundary normal velocity caused by a parameter velocity "deltaP"
    // TODO: update to make per-vertex effectiveVelocity (instead of
    // effective normal velocity)
    SField effectiveNormalVelocity(const SField &deltaP) const {
        SField vn(mesh().numBoundaryElements());
        for (size_t bei = 0; bei < vn.size(); ++bei) {
            vn[bei] = 0;
            for (size_t p = 0; p < deltaP.size(); ++p)
                vn[bei] += deltaP[p] * m_vn_p[p][bei].average();
        }
        return vn;
    }

    const NormalShapeVelocity<_N> &parameterNormalVelocities(size_t p) const { return m_vn_p[p]; }
    const std::vector<VField    > &fluctuationDisplacements()          const { return w_ij; }

    void writeDescription(std::ostream &os) const {
        os << "p:";
        for (size_t i = 0; i < m_params.size(); ++i)
            os << "\t" << m_params[i];
        os << std::endl;

        os << "moduli:\t";
        C.printOrthotropic(os);
        os << "anisotropy:\t" << C.anisotropy() << std::endl;
        os << "JS:\t" << evaluateJS() << std::endl;
        os << "printable:\t" << m_printable << std::endl;

        SField gradP = gradp_JS();
        os << "grad_p(J_S):\t";
        gradP.print(os, "", "", "", "\t");
        os << std::endl << "||grad_p Js||:\t" << gradP.norm() << std::endl;
    }

    VField directionField(const SField &v_n) const {
        size_t numBE = mesh().numBoundaryElements();
        assert(v_n.domainSize() == numBE);
        VField direction(numBE);
        for (size_t be = 0; be < numBE; ++be)
            direction(be) = v_n[be] * mesh().boundaryElement(be)->normal();
        return direction;
    }

    void writeMeshAndFields(const std::string &name) const {
        auto complianceFitGrad = shapeDerivativeJS();
        SField avg_vn(complianceFitGrad.size());
        for (size_t i = 0; i < complianceFitGrad.size(); ++i)
            avg_vn[i] = complianceFitGrad[i].average();
        auto projectedNormalVelocity = effectiveNormalVelocity(gradp_JS());

        if (_N == 2) {
            MeshIO::save(name + ".msh", mesh());
            EdgeFields ef(mesh());
            ef.addField("avg_gradFit", avg_vn);
            ef.addField("avg_gradFit_direction", directionField(avg_vn));

            ef.addField("projectedVn", projectedNormalVelocity);
            ef.addField("projectedVn_direction", directionField(projectedNormalVelocity));
            ef.write(name + ".ef");
        }
        if (_N == 3) {
            std::vector<MeshIO::IOVertex>  bdryVertices;
            std::vector<MeshIO::IOElement> bdryElements;
            const auto &m = mesh();
            for (size_t i = 0; i < m.numBoundaryVertices(); ++i) {
                bdryVertices.emplace_back(m.boundaryVertex(i).volumeVertex().node()->p);
            }
            for (size_t i = 0; i < m.numBoundaryElements(); ++i) {
                auto be = m.boundaryElement(i);
                bdryElements.emplace_back(be.vertex(0).index(),
                                          be.vertex(1).index(),
                                          be.vertex(2).index());
            }

            MSHFieldWriter writer(name + ".msh", bdryVertices, bdryElements);
            for (size_t p = 0; p < m_vn_p.size(); ++p) {
                const auto &vn = m_vn_p[p];;
                SField nvel(m.numBoundaryElements());
                for (size_t bei = 0; bei < m.numBoundaryElements(); ++bei)
                    nvel[bei] = vn[bei].average();
                writer.addField("avg_normal_velocity" + std::to_string(p), nvel, DomainType::PER_ELEMENT);
                writer.addField("avg_normal_velocity_direction" + std::to_string(p), directionField(nvel), DomainType::PER_ELEMENT);
            }
            writer.addField("avg_gradFit", avg_vn, DomainType::PER_ELEMENT);
            writer.addField("avg_gradFit_direction", directionField(avg_vn), DomainType::PER_ELEMENT);

            writer.addField("projectedVn", projectedNormalVelocity, DomainType::PER_ELEMENT);
            writer.addField("projectedVn_direction", directionField(projectedNormalVelocity), DomainType::PER_ELEMENT);
        }
    }

    void writeVolumeMesh(const std::string &name) const {
        MeshIO::save(name, mesh());
    }

    void dumpSimulationMatrix(const std::string &matOut) const {
        m_sim->dumpSystem(matOut);
    }

    // Note, must overwrite inflator's parameter state :(
    template<class _Inflator>
    void writePatternDoFs(const std::string &name, _Inflator &inflator) {
        inflator.writePatternDoFs(name, m_params);
    }

    bool paramsDiffer(size_t nParams, const Real *params) const {
        assert(nParams = m_params.size());
        for (size_t i = 0; i < nParams; ++i)
            if (m_params[i] != params[i])
                return true;
        return false;
    }

    const _Sim &simulator() const { return *m_sim; }
          _Sim &simulator()       { return *m_sim; }

          typename _Sim::Mesh &mesh()       { return m_sim->mesh(); }
    const typename _Sim::Mesh &mesh() const { return m_sim->mesh(); }

    const _ETensor &elasticityTensor() const { return C; }
    const _ETensor &complianceTensor() const { return S; }

    ////////////////////////////////////////////////////////////////////////////
    // Target volume objective implementation
    ////////////////////////////////////////////////////////////////////////////
    void setTargetVolume(Real v) { m_targetVol = v; }

    // Shape derivative of the volume functional:
    // Vol(omega) = int_omega 1 dV  ==> dVol[v] = int_domega 1 * (v . n) dA
    typedef Interpolant<Real, BEGradTensorInterpolant::K, 0> BEConstInterpolant;
    std::vector<BEConstInterpolant> shapeDerivativeVolume() const {
        std::vector<BEConstInterpolant> grad(m_gradS.size());
        for (size_t be = 0; be < m_gradS.size(); ++be)
            grad[be] = 1.0;
        return grad;
    }

    // Evaluate (unweighted) volume-fitting objective weight.
    // ||V - V^*||^2
    Real evaluateJVol() const {
        if (!m_targetVol) throw std::runtime_error("Target volume not set.");
        Real diff = *m_targetVol - mesh().volume();
        return diff * diff;
    }

    // Shape derivative of (unweighted) volume-fitting objective.
    // JVol(omega) = ||Vol - V^*||^2 => dJVol[v] = 2 (Vol - V^*) dVol[v]
    std::vector<BEConstInterpolant> shapeDerivativeJVol() const {
        // TODO: make velocity zero on periodic elements?
        // Shouldn't be necessary for most inflators, though, since they enforce
        // zero parameter velocity on the periodic boundary.
        if (!m_targetVol) throw std::runtime_error("Target volume not set.");
        std::vector<BEConstInterpolant> grad(m_gradS.size());
        Real diff = mesh().volume() - *m_targetVol;
        for (size_t be = 0; be < m_gradS.size(); ++be)
            grad[be] = 2 * diff;
        return grad;
    }

    // Shape derivative of (unweighted) volume-fitting objective.
    // JVol(omega) = ||Vol - V^*||^2 => dJVol[v] = 2 (Vol - V^*) dVol[v]
    SField gradp_JVol() const {
        if (!m_targetVol) throw std::runtime_error("Target volume not set.");

        SField result(m_vn_p.size());
        result.clear();

        // First compute gradient of volume: int v_n dV
        for (size_t p = 0; p < m_vn_p.size(); ++p) {
            for (auto be : mesh().boundaryElements()) {
                if (be->isPeriodic) continue;
                result[p] += m_vn_p[p][be.index()].integrate(be->volume());
            }
        }

        Real diff = mesh().volume() - *m_targetVol;
        result *= (2.0 * diff);
        return result;
    }

    bool printable() const { return m_printable; }

protected:
    std::unique_ptr<_Sim> m_sim;
    _ETensor C, S, m_targetS, m_diffS;
    boost::optional<Real> m_targetVol;
    std::vector<BEGradTensorInterpolant> m_gradS;
    std::vector<_ETensor>                m_gradp_S;
    std::vector<NormalShapeVelocity<_N>> m_vn_p;
    bool m_printable;

    // Fluctuation displacements--only kept if requested in constructor (a subclasses
    // might need them).
    std::vector<VField> w_ij;

    // Requests linear objective/residual estimate for when meshing fails
    std::vector<Real> m_estimateObjectiveWithDeltaP;
    std::vector<Real> m_params;
};

// Use previous iterate if evaluating the same point. Otherwise, attempt to
// inflate the new parameters. Try three times to inflate, and if unsuccessful,
// estimate the point with linear extrapolation.
//
// Because this is a variadic template, it can be used with an iterate class
// whose constructor has arbitrary trailing arguments.
template<class _Iterate, class _Inflator, typename... Args>
std::unique_ptr<_Iterate>
getIterate(std::unique_ptr<_Iterate> oldIterate,
        _Inflator &inflator, size_t nParams, const double *params,
        Args&&... args) {
    if (!oldIterate || oldIterate->paramsDiffer(nParams, params)) {
        std::unique_ptr<_Iterate> newIterate;
        bool success;
        for (size_t i = 0; i < 3; ++i) {
            success = true;
            try {
                newIterate = Future::make_unique<_Iterate>(inflator,
                                nParams, params, std::forward<Args>(args)...);
            }
            catch (std::exception &e) {
                std::cerr << "INFLATOR FAILED: " << e.what() << std::endl;
                success = false;
            }
            if (success) break;
        }
        if (!success) {
            std::cerr << "3 INFLATION ATTEMPTS FAILED." << std::endl;
            if (!oldIterate) throw std::runtime_error("Inflation failure on first iterate");
            newIterate = std::move(oldIterate);
            newIterate->estimatePoint(nParams, params);
        }
        return newIterate;
    }
    else {
        // Old iterate is exact, not an approximation
        oldIterate->disableEstimation();
        return oldIterate;
    }
}

}

#endif /* end of include guard: PATTERNOPTIMIZATIONITERATE_HH */
