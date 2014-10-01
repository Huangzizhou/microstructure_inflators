////////////////////////////////////////////////////////////////////////////////
// PatternOptimization.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Nonlinear least squares-based pattern parameter optimizer.
//      Attempts to fit compliance tensors:
//          1/2 ||S - S^*||^2_F
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/29/2014 14:22:30
////////////////////////////////////////////////////////////////////////////////
#ifndef PATTERNOPTIMIZATION_HH
#define PATTERNOPTIMIZATION_HH

#include <ceres/ceres.h>
#include <glog/logging.h>

#include <cassert>
#include <memory>

#include <WireInflator2D.h>

namespace PatternOptimization {

template<size_t _N>
class Optimizer {
    typedef typename LinearElasticityND<_N>:: template
        HomogenousSimulator<Materials::Constant>       Simulator;
    typedef typename LinearElasticityND<_N>::VField VField;
    typedef ScalarField<Real> SField;
    typedef typename LinearElasticityND<_N>::ETensor _ETensor;
public:
    Optimizer(WireInflator2D &inflator,
            const TessellationParameters &tparams)
        : m_inflator(inflator), m_tparams(tparams) {
    }

    struct Iterate {
        typedef typename LinearElasticityND<_N>:: template
            HomogenousSimulator<Materials::Constant>       Simulator;
        Iterate(WireInflator2D &inflator, const TessellationParameters &tparams,
                size_t nParams, const double *params, const _ETensor &targetS)
            : m_sim(NULL), m_targetS(targetS)
        {
            m_params.resize(nParams);
            CellParameters p_params = inflator.createParameters();
            for (size_t i = 0; i < m_params.size(); ++i) {
                m_params[i] = params[i];
                p_params.parameter(i) = params[i];
                std::cout << params[i] << "\t";
            }
            std::cout << endl;

            if (!inflator.patternGenerator().parametersValid(p_params))
                throw runtime_error("Invalid parameters specified.");
            
            WireInflator2D::OutMeshType inflatedMesh;
            inflator.generatePattern(p_params, tparams, inflatedMesh);
            m_sim = new Simulator(inflatedMesh.elements, inflatedMesh.nodes);

            size_t numBE = m_sim->mesh().numBoundaryElements();
            vn_p.assign(nParams, SField(numBE));
            for (size_t bei = 0; bei < numBE; ++bei) {
                auto be = m_sim->mesh().boundaryElement(bei);
                auto edge = make_pair(be.tip(). volumeVertex().index(),
                                      be.tail().volumeVertex().index());
                const auto &field = inflatedMesh.edge_fields.at(edge);
                assert(field.size() == nParams);
                for (size_t p = 0; p < nParams; ++p)
                    vn_p[p][bei] = field[p];
            }

            std::vector<VField> w_ij;
            PeriodicHomogenization::solveCellProblems(w_ij, *m_sim);
            C = PeriodicHomogenization::homogenizedElasticityTensor(w_ij, *m_sim);
            S = C.inverse();
            std::vector<_ETensor> gradEh =
                PeriodicHomogenization::homogenizedTensorGradient(w_ij, *m_sim);
            for (const auto &G : gradEh)
                gradS.push_back(-S.doubleDoubleContract(G));
        }

        // Evaluate compliance frobenius norm objective.
        Real evaluateJS() const {
            auto diff = S - m_targetS;
            return 0.5 * diff.quadrupleContract(diff);
        }

        ////////////////////////////////////////////////////////////////////////
        /*! Computes grad(1/2 sum_ijkl (S_ijkl - target_ijlk|)^2) =
        //      (S_ijkl - target_ijlk) * grad(S_ikjl))
        //  @param[in]  target  S^* (target compliance tensor)
        //  @return     SField  Per-boundary-edge scalar field giving steepest
        //                      ascent normal velocity perturbation for JS
        *///////////////////////////////////////////////////////////////////////
        SField shapeDerivativeJS() const {
            _ETensor diff = S - m_targetS;
            SField v_n(gradS.size());

            for (size_t be = 0; be < gradS.size(); ++be)
                v_n[be] = diff.quadrupleContract(gradS[be]);

            return v_n;
        }

        SField gradp_JS(const SField &gradJS) const {
            SField result(vn_p.size());
            result.clear();
            for (size_t p = 0; p < vn_p.size(); ++p) {
                for (size_t bei = 0; bei < m_sim->mesh().numBoundaryElements(); ++bei) {
                    auto be = m_sim->mesh().boundaryElement(bei);
                    result[p] += vn_p[p][bei] * gradJS[bei] * be->area();
                }
            }
            return result;
        }
        SField gradp_JS() const { return gradp_JS(shapeDerivativeJS()); }

        // The (i, j)th residual (j >= i) for the nonlinear least squares (a
        // single term of the Frobenius distance). The terms are weighted so
        // that the squared norm of the residual vector corresponds to the
        // Frobenius norm of the rank 4 tensor difference S - S^*.
        Real residual(size_t i, size_t j) const {
            assert(j >= i);
            Real weight = 1.0;
            if (j != i)  weight *= sqrt(2); // Account for lower triangle
            if (i >= _N) weight *= sqrt(2); // Left shear doubler
            if (j >= _N) weight *= sqrt(2); // Right shear doubler
            return weight * (S.D(i, j) - m_targetS.D(i, j));
        }

        // Derivative of residual(i, j) wrt parameter p:
        // d/dp (S_ijkl - target_ijkl) = d/dp S_ijkl = <gradS_ijkl, vn_p>
        Real jacobian(size_t i, size_t j, size_t p) const {
            assert(j >= i);
            Real result = 0;
            for (size_t bei = 0; bei < m_sim->mesh().numBoundaryElements(); ++bei) {
                auto be = m_sim->mesh().boundaryElement(bei);
                result += vn_p[p][bei] * gradS[bei].D(i, j) * be->area();
            }

            Real weight = 1.0;
            if (j != i)  weight *= sqrt(2); // Account for lower triangle
            if (i >= _N) weight *= sqrt(2); // Left shear doubler
            if (j >= _N) weight *= sqrt(2); // Right shear doubler
            return weight * result;
        }

        // Boundary normal velocity caused by a parameter velocity "deltaP"
        SField effectiveNormalVelocity(const SField &deltaP) const {
            SField vn(m_sim->mesh().numBoundaryElements());
            for (size_t bei = 0; bei < vn.size(); ++bei) {
                vn[bei] = 0;
                for (size_t p = 0; p < deltaP.size(); ++p)
                    vn[bei] += deltaP[p] * vn_p[p][bei];
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
            os << "residual JS:\t";

            SField gradP = gradp_JS();
            os << "grad_p(J_S):\t";
            gradP.print(os, "", "", "", "\t");
            os << std::endl << "||grad_p)||:\t" << gradP.norm() << std::endl;
        }

        VField directionField(const SField &v_n) const {
            size_t numBE = m_sim->mesh().numBoundaryElements();
            VField direction(numBE);
            for (size_t be = 0; be < numBE; ++be)
                direction(be) = v_n[be] * m_sim->mesh().boundaryElement(be)->normal();
            return direction;
        }

        void writeMeshAndFields(const std::string &name) const {
            MeshIO::save(name + ".msh", m_sim->mesh());
            EdgeFields ef(m_sim->mesh());
            auto complianceFitGrad = shapeDerivativeJS();
            ef.addField("gradFit", complianceFitGrad);
            ef.addField("gradFit direction", directionField(complianceFitGrad));

            auto projectedNormalVelocity =
                effectiveNormalVelocity(gradp_JS(complianceFitGrad));
            ef.addField("projectedVn", projectedNormalVelocity);
            ef.addField("projectedVn direction", directionField(projectedNormalVelocity));
            ef.write(name + ".ef");
        }

        bool paramsDiffer(size_t nParams, Real *params) const {
            for (size_t i = 0; i < nParams; ++i)
                if (m_params[i] != params[i])
                    return true;
            return false;
        }

    private:
        Simulator *m_sim;
        _ETensor C, S, m_targetS;
        std::vector<_ETensor> gradS;

        // Parameter normal velocity fields
        std::vector<SField> vn_p;

        std::vector<Real> m_params;
    };

    struct TensorFitCost : public ceres::CostFunction {
        typedef ceres::CostFunction Base;
        TensorFitCost(WireInflator2D &inflator,
                const TessellationParameters &tparams, const _ETensor &targetS)
            : m_inflator(inflator), m_tparams(tparams), m_targetS(targetS) {
            Base::set_num_residuals((_N == 2) ? 6 : 21);
            // We put all the pattern parameters in a single parameter block.
            Base::mutable_parameter_block_sizes()->assign(1,
                    m_inflator.patternGenerator().numberOfParameters());
        }

        virtual bool Evaluate(double const * const *parameters,
                double *residuals, double **jacobians) const {
            size_t nParams = parameter_block_sizes()[0];
            Iterate it(m_inflator, m_tparams, nParams,
                       parameters[0], m_targetS);

            size_t r = 0;
            for (size_t i = 0; i < flatLen(_N); ++i)
                for (size_t j = i; j < flatLen(_N); ++j)
                    residuals[r++] = it.residual(i, j);

            if (jacobians == NULL) return true;

            r = 0;
            for (size_t i = 0; i < flatLen(_N); ++i) {
                for (size_t j = i; j < flatLen(_N); ++j) {
                    for (size_t p = 0; p < nParams; ++p)
                        jacobians[0][r * nParams + p] = it.jacobian(i, j, p);
                    ++r;
                }
            }

            return true;
        }

        WireInflator2D &m_inflator;
        TessellationParameters m_tparams;
        _ETensor m_targetS;
    };

    void optimize(SField &params, const _ETensor &targetS,
                  const string outName) {
        ceres::CostFunction *fitCost =
                new TensorFitCost(m_inflator, m_tparams, targetS);
        ceres::Problem problem;
        problem.AddResidualBlock(fitCost, NULL, params.data());

        for (size_t p = 0; p < params.domainSize(); ++p) {
            auto range = m_inflator.patternGenerator().getParameterRange(p);
            std::cout << "setting bound on variable " << p << ": " << range.first << ", " << range.second << std::endl;
            problem.SetParameterLowerBound(params.data(), p, range.first);
            problem.SetParameterUpperBound(params.data(), p, range.second);
        }

        ceres::Solver::Options options;
        // options.minimizer_type = ceres::LINE_SEARCH;
        // options.line_search_direction_type = ceres::BFGS;
        // options.trust_region_strategy_type = ceres::DOGLEG;
        // options.dogleg_type = ceres::SUBSPACE_DOGLEG;
        options.use_nonmonotonic_steps = true;
        options.minimizer_progress_to_stdout = true;
        options.initial_trust_region_radius = 0.01;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        cout << summary.BriefReport() << "\n";
    }

private:
    WireInflator2D &m_inflator;
    TessellationParameters m_tparams;
    std::vector<Real> m_patternParams;
};

}

#endif /* end of include guard: PATTERNOPTIMIZATION_HH */
