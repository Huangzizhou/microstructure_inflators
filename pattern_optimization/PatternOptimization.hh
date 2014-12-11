////////////////////////////////////////////////////////////////////////////////
// PatternOptimization.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Nonlinear least squares-based pattern parameter optimizer.
//      Attempts to fit compliance tensors:
//          1/2 ||S - S^*||^2_F
//      Note: should only be used with a homogenous material-backed simulator.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  09/29/2014 14:22:30
////////////////////////////////////////////////////////////////////////////////
#ifndef PATTERNOPTIMIZATION_HH
#define PATTERNOPTIMIZATION_HH

#include <ceres/ceres.h>
#include <glog/logging.h>
#include <dlib/optimization.h>

#include <cassert>
#include <memory>

#include <WireInflator2D.h>
#include <EdgeFields.hh>

namespace PatternOptimization {

template<class _Sim>
class Optimizer {
    typedef typename _Sim::VField VField;
    typedef ScalarField<Real> SField;
    typedef typename _Sim::ETensor _ETensor;
    typedef dlib::matrix<double,0,1> dlib_vector;
    static constexpr size_t _N = _Sim::N;
public:
    Optimizer(WireInflator2D &inflator,
            const TessellationParameters &tparams, std::vector<Real> radiusBounds,
            std::vector<Real> translationBounds)
        : m_inflator(inflator), m_tparams(tparams), m_radiusBounds(radiusBounds),
          m_transBounds(translationBounds) { }

    template<class _Vector>
    void getParameterBounds(_Vector &lowerBounds, _Vector &upperBounds) {
        auto ops = m_inflator.patternGenerator().getParameterOperations();
        for (size_t p = 0; p < ops.size(); ++p) {
            switch (ops.at(p).type) {
                case ParameterOperation::Radius:
                    lowerBounds(p) = m_radiusBounds.at(0);
                    upperBounds(p) = m_radiusBounds.at(1);
                    break;
                case ParameterOperation::Translation:
                    lowerBounds(p) = m_transBounds.at(0);
                    upperBounds(p) = m_transBounds.at(1);
                    break;
                default: assert(false);
            }
        }
    }

    struct Iterate {
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

        Iterate(WireInflator2D &inflator, const TessellationParameters &tparams,
                size_t nParams, const double *params, const _ETensor &targetS)
            : m_targetS(targetS)
        {
            m_params.resize(nParams);
            CellParameters p_params = inflator.createParameters();
            for (size_t i = 0; i < m_params.size(); ++i) {
                m_params[i] = params[i];
                p_params.parameter(i) = params[i];
            }

            if (!inflator.patternGenerator().parametersValid(p_params))
                throw runtime_error("Invalid parameters specified.");
            
            WireInflator2D::OutMeshType inflatedMesh;
            inflator.generatePattern(p_params, tparams, inflatedMesh);
            m_sim = std::make_shared<_Sim>(inflatedMesh.elements,
                                           inflatedMesh.nodes);

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
        //  @return Per-boundary-edge piecewise 2 * (deg - 1) scalar field giving
        //          steepest ascent normal velocity perturbation for JS
        *///////////////////////////////////////////////////////////////////////
        std::vector<BEGradInterpolant> shapeDerivativeJS() const {
            _ETensor diff = S - m_targetS;
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
            SField result(vn_p.size());
            result.clear();
            for (size_t p = 0; p < vn_p.size(); ++p) {
                for (size_t bei = 0; bei < m_sim->mesh().numBoundaryElements(); ++bei) {
                    auto be = m_sim->mesh().boundaryElement(bei);
                    Real vn = vn_p[p][bei]; // TODO: make this a linear interpolant.
                    const auto &grad = gradJS[bei];
                    result[p] +=
                        Quadrature<_Sim::K - 1, 1 + BEGradTensorInterpolant::Deg>::
                        integrate([&] (const VectorND<be.numVertices()> &pt) {
                            return vn * grad(pt);
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
                Real vn = vn_p[p][bei]; // TODO: make this a linear interpolant.
                const auto &grad = m_gradS[bei];
                result +=
                    Quadrature<_Sim::K - 1, 1 + BEGradTensorInterpolant::Deg>::
                    integrate([&] (const VectorND<be.numVertices()> &pt) {
                        return vn * grad(pt).D(ij, jk);
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
            MeshIO::save(name + ".msh", m_sim->mesh());
            EdgeFields ef(m_sim->mesh());
            auto complianceFitGrad = shapeDerivativeJS();
            SField avg_vn(complianceFitGrad.size());
            for (size_t i = 0; i < complianceFitGrad.size(); ++i)
                avg_vn[i] = complianceFitGrad[i].average();
            ef.addField("(averaged) gradFit", avg_vn);
            ef.addField("(averaged) gradFit direction", directionField(avg_vn));

            auto projectedNormalVelocity =
                effectiveNormalVelocity(gradp_JS(complianceFitGrad));
            ef.addField("projectedVn", projectedNormalVelocity);
            ef.addField("projectedVn direction", directionField(projectedNormalVelocity));
            ef.write(name + ".ef");
        }

        bool paramsDiffer(size_t nParams, const Real *params) const {
            assert(nParams = m_params.size());
            for (size_t i = 0; i < nParams; ++i)
                if (m_params[i] != params[i])
                    return true;
            return false;
        }

    private:
        std::shared_ptr<_Sim> m_sim;
        _ETensor C, S, m_targetS;
        std::vector<BEGradTensorInterpolant> m_gradS;

        // Parameter normal velocity fields
        std::vector<SField> vn_p;

        std::vector<Real> m_params;
    };

    // Forward declaration so friend-ing can happen
    class IterationCallback;

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
            // Ceres seems to like re-evaluating at the same point, so we detect
            // this to avoid re-solving.
            if (!m_iterate || m_iterate->paramsDiffer(nParams, parameters[0])) {
                m_iterate = std::make_shared<Iterate>(m_inflator, m_tparams,
                        nParams, parameters[0], m_targetS);
            }

            size_t r = 0;
            for (size_t i = 0; i < flatLen(_N); ++i)
                for (size_t j = i; j < flatLen(_N); ++j)
                    residuals[r++] = m_iterate->residual(i, j);

            if (jacobians == NULL) return true;

            r = 0;
            for (size_t i = 0; i < flatLen(_N); ++i) {
                for (size_t j = i; j < flatLen(_N); ++j) {
                    for (size_t p = 0; p < nParams; ++p)
                        jacobians[0][r * nParams + p] = m_iterate->jacobian(i, j, p);
                    ++r;
                }
            }

            return true;
        }

        std::shared_ptr<const Iterate> currentIterate() const { return m_iterate; }

        virtual ~TensorFitCost() { }

    private:
        WireInflator2D &m_inflator;
        TessellationParameters m_tparams;
        _ETensor m_targetS;
        // Ceres requires Evaluate to be constant, so this caching pointer must
        // be made mutable.
        mutable std::shared_ptr<Iterate> m_iterate;

        friend class IterationCallback;
    };

    class IterationCallback : public ceres::IterationCallback {
    public:
        IterationCallback(TensorFitCost &evalulator, SField &params,
                const std::string &outPath)
            : m_evaluator(evalulator), m_params(params), m_outPath(outPath),
            m_iter(0) {}
        ceres::CallbackReturnType operator()(const ceres::IterationSummary &sum)
        {
            size_t nParams = m_params.domainSize();
            auto curr = m_evaluator.currentIterate();
            if (curr->paramsDiffer(nParams, &m_params[0])) {
                curr = std::make_shared<Iterate>(m_evaluator.m_inflator,
                        m_evaluator.m_tparams, nParams, &m_params[0],
                        m_evaluator.m_targetS);
            }
            curr->writeDescription(std::cout);
            std::cout << std::endl;

            if (m_outPath != "")
                curr->writeMeshAndFields(std::to_string(m_iter) + "_" + m_outPath);

            ++m_iter;
            return ceres::SOLVER_CONTINUE;
        }

        virtual ~IterationCallback() { }
    private:
        TensorFitCost &m_evaluator;
        SField &m_params;
        std::string m_outPath;
        size_t m_iter;
    };

    void optimize_lm(SField &params, const _ETensor &targetS,
                  const string &outPath) {
        TensorFitCost *fitCost = new TensorFitCost(m_inflator, m_tparams,
                                                   targetS);
        ceres::Problem problem;
        problem.AddResidualBlock(fitCost, NULL, params.data());

        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        for (size_t p = 0; p < params.domainSize(); ++p) {
            problem.SetParameterLowerBound(params.data(), p, lowerBounds[p]);
            problem.SetParameterUpperBound(params.data(), p, upperBounds[p]);
        }

        ceres::Solver::Options options;
        options.update_state_every_iteration = true;
        IterationCallback cb(*fitCost, params, outPath);
        options.callbacks.push_back(&cb);
        // options.minimizer_type = ceres::LINE_SEARCH;
        // options.line_search_direction_type = ceres::BFGS;
        // options.trust_region_strategy_type = ceres::DOGLEG;
        // options.dogleg_type = ceres::SUBSPACE_DOGLEG;
        // options.use_nonmonotonic_steps = true;
        // options.minimizer_progress_to_stdout = true;
        // options.initial_trust_region_radius = 0.01;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << "\n";
    }

    void optimize_gd(SField &params, const _ETensor &targetS,
            size_t niters, double alpha, const string &outName) {
        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        // Gradient Descent Version
        for (size_t i = 0; i < niters; ++i) {
            Iterate iterate(m_inflator, m_tparams, params.size(), params.data(), targetS);
            SField gradP = iterate.gradp_JS();
            iterate.writeDescription(std::cout);
            std::cout << std::endl;

            params -= gradP * alpha;

            // Apply bound constraints
            params.minRelax(upperBounds);
            params.maxRelax(lowerBounds);

            if (outName != "")
                iterate.writeMeshAndFields(std::to_string(i) + "_" + outName);
        }
    }

    // Evaluates the objective by inflating the wire mesh and homogenizing.
    // The iterate stored internally also knows how to evaluate the gradient
    // efficiently, so our GradientEvaluator below just accesses it.
    struct DLibObjectiveEvaluator {
        DLibObjectiveEvaluator(WireInflator2D &inflator,
                const TessellationParameters &tparams, const _ETensor &targetS)
            : m_iterate(NULL), m_inflator(inflator), m_tparams(tparams),
              m_targetS(targetS) {
            nParams = m_inflator.patternGenerator().numberOfParameters();
        };

        double operator()(const dlib_vector &x) const {
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);

            m_iterate = std::make_shared<Iterate>(m_inflator, m_tparams,
                    nParams, &x_vec[0], m_targetS);
            return m_iterate->evaluateJS();
        }

        const Iterate &currentIterate() const {
            assert(m_iterate); return *m_iterate;
        };

        size_t nParams;
    private:
        // Iterate is mutable so that operator() can be const as dlib requires
        mutable std::shared_ptr<Iterate> m_iterate;
        WireInflator2D &m_inflator;
        TessellationParameters m_tparams;
        _ETensor m_targetS;
    };
    
    // Extracts gradient from the iterate constructed by DLibObjectiveEvaluator.
    struct DLibGradientEvaluator {
        DLibGradientEvaluator(const DLibObjectiveEvaluator &obj) : m_obj(obj) { }

        dlib_vector operator()(const dlib_vector &x) const {
            size_t nParams = m_obj.nParams;
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);
            if (m_obj.currentIterate().paramsDiffer(nParams, &x_vec[0]))
                throw std::runtime_error("Objective must be evaluated first");
            SField gradp(m_obj.currentIterate().gradp_JS());

            dlib_vector result(nParams);
            for (size_t p = 0; p < nParams; ++p)
                result(p) = gradp[p];
            return result;
        }
    private:
        const DLibObjectiveEvaluator &m_obj;
    };

    // Hack to get notified at the end of each iteration: subclass the stop
    // strategy.
    class ReportingStopStrategy : public dlib::objective_delta_stop_strategy {
        typedef dlib::objective_delta_stop_strategy Base;
    public:
        ReportingStopStrategy(double min_delta, unsigned long max_iter,
                              DLibObjectiveEvaluator &obj, const std::string &outPath)
            : Base(min_delta, max_iter), m_obj(obj), m_outPath(outPath), m_iter(0) { }

        template <typename T>
        bool should_continue_search(const T& x, const double funct_value,
            const T& funct_derivative) {
            m_obj.currentIterate().writeDescription(std::cout);
            cout << endl;
            if (m_outPath != "")
                m_obj.currentIterate().writeMeshAndFields(std::to_string(m_iter)
                        + "_" + m_outPath);
            ++m_iter;
            return Base::should_continue_search(x, funct_value, funct_derivative);
        }

    private:
        const DLibObjectiveEvaluator &m_obj;
        std::string m_outPath;
        size_t m_iter;
    };

    // If max_size = 0, plain bfgs is used
    // otherwise l-bfgs is used.
    void optimize_bfgs(SField &params, const _ETensor &targetS, size_t niters,
                       const string &outPath, size_t max_size = 0) {
        DLibObjectiveEvaluator obj(m_inflator, m_tparams, targetS);
        DLibGradientEvaluator grad(obj);

        size_t nParams = m_inflator.patternGenerator().numberOfParameters();
        // convert initial parameter vector
        dlib_vector optParams(nParams);
        for (size_t p = 0; p < nParams; ++p)
            optParams(p) = params[p];

        dlib_vector lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);

        if (max_size == 0)
            dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
                    ReportingStopStrategy(1e-16, niters, obj, outPath),
                    obj, grad, optParams, lowerBounds, upperBounds);
        else
            dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(max_size),
                    ReportingStopStrategy(1e-16, niters, obj, outPath),
                    obj, grad, optParams, lowerBounds, upperBounds);

        // convert solution
        for (size_t p = 0; p < nParams; ++p)
            params[p] = optParams(p);
    }

private:
    WireInflator2D &m_inflator;
    TessellationParameters m_tparams;
    std::vector<Real> m_patternParams, m_radiusBounds, m_transBounds;
};

}

#endif /* end of include guard: PATTERNOPTIMIZATION_HH */
