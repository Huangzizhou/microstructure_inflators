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

#include "Inflator.hh"
#include <EdgeFields.hh>
#include <MSHFieldWriter.hh>

#include <cassert>
#include <memory>

#include <MeshIO.hh>

namespace PatternOptimization {
template<class _Sim>
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
    // represent the per-boundary-element value of the shape derivative a
    // scalar function.
    typedef Interpolant<Real, BEGradTensorInterpolant::K,
                    BEGradTensorInterpolant::Deg>   BEGradInterpolant;

    Iterate(Inflator<_N> &inflator, size_t nParams, const double *params,
            const _ETensor &targetS)
        : m_targetS(targetS)
    {
        m_params.resize(nParams);
        for (size_t i = 0; i < m_params.size(); ++i)
            m_params[i] = params[i];
        BENCHMARK_START_TIMER("Inflate");
        inflator.inflate(m_params);
        BENCHMARK_STOP_TIMER("Inflate");

        BENCHMARK_START_TIMER_SECTION("Eval");
        m_sim = std::make_shared<_Sim>(inflator.elements(),
                                       inflator.vertices());
        m_vn_p = inflator.computeShapeNormalVelocities(m_sim->mesh());

        std::vector<VField> w_ij;
        PeriodicHomogenization::solveCellProblems(w_ij, *m_sim);
        C = PeriodicHomogenization::homogenizedElasticityTensor(w_ij, *m_sim);
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

        BENCHMARK_STOP_TIMER_SECTION("Eval");
    }

    // Evaluate compliance frobenius norm objective.
    Real evaluateJS() const {
        auto diff = S - m_targetS;
        // TODO: FIGURE OUT WHAT TO DO WITH THIS HACK.
        // Ignore "non-orthotropic part" (this part is erroneous)
        for (size_t j = _N; j < flatLen(_N); ++j) {
            for (size_t i = 0; i < j; ++i) {
                diff.D(i, j) = 0.0;
            }
        }
        return 0.5 * diff.quadrupleContract(diff);
    }

    ////////////////////////////////////////////////////////////////////////
    /*! Computes grad(1/2 sum_ijkl (S_ijkl - target_ijlk|)^2) =
    //      (S_ijkl - target_ijlk) * grad(S_ikjl))
    //  @param[in]  target  S^* (target compliance tensor)
    //  @return Per-boundary-edge piecewise 2 * (deg - 1) scalar field giving
    //          steepest ascent normal velocity perturbation for JS
    *///////////////////////////////////////////////////////////////////////
    std::vector<BEGradInterpolant> shapeDerivativeJS() const {
        _ETensor diff = S - m_targetS;
        // TODO: FIGURE OUT WHAT TO DO WITH THIS HACK.
        // Ignore "non-orthotropic part" (this part is erroneous)
        for (size_t j = _N; j < flatLen(_N); ++j) {
            for (size_t i = 0; i < j; ++i) {
                diff.D(i, j) = 0.0;
            }
        }
        std::vector<BEGradInterpolant> grad(m_gradS.size());

        for (size_t be = 0; be < m_gradS.size(); ++be) {
            // Compute each nodal value of the interpolant.
            const auto &GS = m_gradS[be];
                  auto &g = grad[be];
            for (size_t n = 0; n < GS.size(); ++n)
                g[n] = diff.quadrupleContract(GS[n]);
        }

        return grad;
    }

    SField gradp_JS(const std::vector<BEGradInterpolant> &gradJS) const {
        SField result(m_params.size());
        result.clear();
        for (size_t p = 0; p < m_params.size(); ++p) {
            for (size_t bei = 0; bei < m_sim->mesh().numBoundaryElements(); ++bei) {
                auto be = m_sim->mesh().boundaryElement(bei);
                const auto &vn = m_vn_p[p][bei]; // parameter normal shape velocity interpolant
                const auto &grad = gradJS[bei];
                result[p] +=
                    Quadrature<_Sim::K - 1, 1 + BEGradTensorInterpolant::Deg>::
                    integrate([&] (const VectorND<be.numVertices()> &pt) {
                        return vn(pt) * grad(pt);
                    }, be->volume());
            }
        }
        return result;
    }

    SField gradp_JS() const { return gradp_JS(shapeDerivativeJS()); }

    // The (ij, jk)th residual (jk >= ij) for the nonlinear least squares (a
    // single term of the Frobenius distance). The terms are weighted so
    // that the squared norm of the residual vector corresponds to the
    // Frobenius norm of the rank 4 tensor difference S - S^*.
    Real residual(size_t ij, size_t jk) const {
        assert(jk >= ij);
        Real weight = 1.0;
        if (jk != ij) weight *= sqrt(2); // Account for lower triangle
        if (ij >= _N) weight *= sqrt(2); // Left shear doubler
        if (jk >= _N) weight *= sqrt(2); // Right shear doubler
        return weight * (S.D(ij, jk) - m_targetS.D(ij, jk));
    }

    // Derivative of residual(ij, jk) wrt parameter p:
    // d/dp (S_ijkl - target_ijkl) = d/dp S_ijkl = <gradS_ijkl, vn_p>
    Real jacobian(size_t ij, size_t jk, size_t p) const {
        assert(jk >= ij);
        Real result = 0;
        for (size_t bei = 0; bei < m_sim->mesh().numBoundaryElements(); ++bei) {
            auto be = m_sim->mesh().boundaryElement(bei);
            const auto &vn = m_vn_p[p][bei]; // parameter normal shape velocity interpolant
            const auto &grad = m_gradS[bei];
            result +=
                Quadrature<_Sim::K - 1, 1 + BEGradTensorInterpolant::Deg>::
                integrate([&] (const VectorND<be.numVertices()> &pt) {
                    return vn(pt) * grad(pt).D(ij, jk);
                }, be->volume());
        }

        Real weight = 1.0;
        if (jk != ij) weight *= sqrt(2); // Account for lower triangle
        if (ij >= _N) weight *= sqrt(2); // Left shear doubler
        if (jk >= _N) weight *= sqrt(2); // Right shear doubler
        return weight * result;
    }

    // Boundary normal velocity caused by a parameter velocity "deltaP"
    // TODO: update to make per-vertex effectiveVelocity (instead of
    // effective normal velocity)
    SField effectiveNormalVelocity(const SField &deltaP) const {
        SField vn(m_sim->mesh().numBoundaryElements());
        for (size_t bei = 0; bei < vn.size(); ++bei) {
            vn[bei] = 0;
            for (size_t p = 0; p < deltaP.size(); ++p)
                vn[bei] += deltaP[p] * m_vn_p[p][bei].average();
        }
        return vn;
    }

    void writeDescription(std::ostream &os) const {
        os << "p:";
        for (size_t i = 0; i < m_params.size(); ++i)
            os << "\t" << m_params[i];
        os << std::endl;

        os << "moduli:\t";
        C.printOrthotropic(os);
        os << "JS:\t" << evaluateJS() << std::endl;

        SField gradP = gradp_JS();
        os << "grad_p(J_S):\t";
        gradP.print(os, "", "", "", "\t");
        os << std::endl << "||grad_p||:\t" << gradP.norm() << std::endl;
    }

    VField directionField(const SField &v_n) const {
        size_t numBE = m_sim->mesh().numBoundaryElements();
        assert(v_n.domainSize() == numBE);
        VField direction(numBE);
        for (size_t be = 0; be < numBE; ++be)
            direction(be) = v_n[be] * m_sim->mesh().boundaryElement(be)->normal();
        return direction;
    }

    void writeMeshAndFields(const std::string &name) const {
        auto complianceFitGrad = shapeDerivativeJS();
        SField avg_vn(complianceFitGrad.size());
        for (size_t i = 0; i < complianceFitGrad.size(); ++i)
            avg_vn[i] = complianceFitGrad[i].average();
        auto projectedNormalVelocity =
            effectiveNormalVelocity(gradp_JS(complianceFitGrad));

        if (_N == 2) {
            MeshIO::save(name + ".msh", m_sim->mesh());
            EdgeFields ef(m_sim->mesh());
            ef.addField("(averaged) gradFit", avg_vn);
            ef.addField("(averaged) gradFit direction", directionField(avg_vn));

            ef.addField("projectedVn", projectedNormalVelocity);
            ef.addField("projectedVn direction", directionField(projectedNormalVelocity));
            ef.write(name + ".ef");
        }
        if (_N == 3) {
            std::vector<MeshIO::IOVertex>  bdryVertices;
            std::vector<MeshIO::IOElement> bdryElements;
            const auto &mesh = m_sim->mesh();
            for (size_t i = 0; i < mesh.numBoundaryVertices(); ++i) {
                bdryVertices.emplace_back(mesh.boundaryVertex(i).volumeVertex().node()->p);
            }
            for (size_t i = 0; i < mesh.numBoundaryElements(); ++i) {
                auto be = mesh.boundaryElement(i);
                bdryElements.emplace_back(be.vertex(0).index(),
                                          be.vertex(1).index(),
                                          be.vertex(2).index());
            }

            MSHFieldWriter writer(name + ".msh", bdryVertices, bdryElements);
            for (size_t p = 0; p < m_vn_p.size(); ++p) {
                const auto &vn = m_vn_p[p];;
                SField nvel(mesh.numBoundaryElements());
                for (size_t bei = 0; bei < mesh.numBoundaryElements(); ++bei)
                    nvel[bei] = vn[bei].average();
                writer.addField("(averaged) normal velocity " + std::to_string(p), nvel, MSHFieldWriter::PER_ELEMENT);
                writer.addField("(averaged) normal velocity direction" + std::to_string(p), directionField(nvel), MSHFieldWriter::PER_ELEMENT);
            }
            writer.addField("(averaged) gradFit", avg_vn, MSHFieldWriter::PER_ELEMENT);
            writer.addField("(averaged) gradFit direction", directionField(avg_vn), MSHFieldWriter::PER_ELEMENT);

            writer.addField("projectedVn", projectedNormalVelocity, MSHFieldWriter::PER_ELEMENT);
            writer.addField("projectedVn direction", directionField(projectedNormalVelocity), MSHFieldWriter::PER_ELEMENT);
        }
    }

    void writeVolumeMesh(const std::string &name) const {
        MeshIO::save(name, m_sim->mesh());
    }

    void dumpSimulationMatrix(const std::string &matOut) const {
        m_sim->dumpSystem(matOut);
    }

    bool paramsDiffer(size_t nParams, const Real *params) const {
        assert(nParams = m_params.size());
        for (size_t i = 0; i < nParams; ++i)
            if (m_params[i] != params[i])
                return true;
        return false;
    }

    _Sim &simulator() { return *m_sim; }
    const _ETensor &elasticityTensor() const { return C; }
    const _ETensor &complianceTensor() const { return S; }

private:
    std::shared_ptr<_Sim> m_sim;
    _ETensor C, S, m_targetS;
    std::vector<BEGradTensorInterpolant> m_gradS;
    std::vector<typename Inflator<_N>::NormalShapeVelocity> m_vn_p;

    std::vector<Real> m_params;
};

}

#endif /* end of include guard: PATTERNOPTIMIZATIONITERATE_HH */