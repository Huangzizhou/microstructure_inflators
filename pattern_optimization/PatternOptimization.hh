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
#include <dlib/optimization.h>

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
    typedef dlib::matrix<double,0,1> dlib_vector;
public:
    Optimizer(WireInflator2D &inflator,
            const TessellationParameters &tparams)
        : m_inflator(inflator), m_tparams(tparams) {
    }

    template<class _Vector>
    void getParameterBounds(_Vector &lowerBounds, _Vector &upperBounds) {
        auto ops = m_inflator.patternGenerator().getParameterOperations();
        for (size_t p = 0; p < ops.size(); ++p) {
            switch (ops.at(p).type) {
                case ParameterOperation::Radius:
                    lowerBounds(p) = 0.05;
                    upperBounds(p) = 0.1;
                    break;
                case ParameterOperation::Translation:
                    lowerBounds(p) = -0.15;
                    upperBounds(p) =  0.15;
                    break;
                default: assert(false);
            }
        }
    }

    struct Iterate {
        typedef typename LinearElasticityND<_N>:: template
            HomogenousSimulator<Materials::Constant>       Simulator;
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
            m_sim = std::make_shared<Simulator>(inflatedMesh.elements,
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

        bool paramsDiffer(size_t nParams, const Real *params) const {
            assert(nParams = m_params.size());
            for (size_t i = 0; i < nParams; ++i)
                if (m_params[i] != params[i])
                    return true;
            return false;
        }

    private:
        std::shared_ptr<Simulator> m_sim;
        _ETensor C, S, m_targetS;
        std::vector<_ETensor> gradS;

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
        IterationCallback(TensorFitCost &evalulator, SField &params)
            : m_evaluator(evalulator), m_params(params) {}
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
            return ceres::SOLVER_CONTINUE;
        }

        virtual ~IterationCallback() { }
    private:
        TensorFitCost &m_evaluator;
        SField &m_params;
    };

    void optimize_lm(SField &params, const _ETensor &targetS,
                  const string outName) {
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
        IterationCallback cb(*fitCost, params);
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
            size_t niters, double alpha, const string outName) {
        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        // Gradient Descent Version
        for (size_t i = 0; i < niters; ++i) {
            typename Optimizer<_N>::Iterate iterate(m_inflator, m_tparams,
                    params.size(), params.data(), targetS);
            SField gradP = iterate.gradp_JS();
            iterate.writeDescription(std::cout);
            std::cout << std::endl;

            params -= gradP * alpha;

            // Apply bound constraints
            params.minRelax(upperBounds);
            params.maxRelax(lowerBounds);

            if (outName != "")
                iterate.writeMeshAndFields(to_string(i) + "_" + outName);
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

    // If max_size = 0, plain bfgs is used
    // otherwise l-bfgs is used.
    void optimize_bfgs(SField &params, const _ETensor &targetS,
                       const string outName, size_t max_size = 0) {
        DLibObjectiveEvaluator obj(m_inflator, m_tparams, targetS);
        DLibGradientEvaluator grad(obj);

        size_t nParams = m_inflator.patternGenerator().numberOfParameters();
        // convert initial parameter vector
        dlib_vector optParams(nParams);
        for (size_t p = 0; p < nParams; ++p)
            optParams(p) = params[p];

        dlib_vector lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);

        // TODO: subclass dlib::objective_delta_stop_strategy to output
        // iterates...
        if (max_size == 0)
            dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
                    dlib::objective_delta_stop_strategy(0).be_verbose(),
                    obj, grad, optParams, lowerBounds, upperBounds);
        else
            dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(max_size),
                    dlib::objective_delta_stop_strategy(0).be_verbose(),
                    obj, grad, optParams, lowerBounds, upperBounds);
        // convert solution
        for (size_t p = 0; p < nParams; ++p)
            params[p] = optParams(p);
    }

private:
    WireInflator2D &m_inflator;
    TessellationParameters m_tparams;
    std::vector<Real> m_patternParams;
};

}

#endif /* end of include guard: PATTERNOPTIMIZATION_HH */
