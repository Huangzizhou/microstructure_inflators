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

#include "Inflator.hh"
#include "PatternOptimizationIterate.hh"
#include <EdgeFields.hh>
#include <MSHFieldWriter.hh>

namespace PatternOptimization {

template<class _Sim>
class Optimizer {
    typedef typename _Sim::VField VField;
    typedef ScalarField<Real> SField;
    typedef typename _Sim::ETensor _ETensor;
    typedef dlib::matrix<double,0,1> dlib_vector;
    static constexpr size_t _N = _Sim::N;
public:
    typedef ::PatternOptimization::Iterate<_Sim> Iterate;
    Optimizer(Inflator<_N> &inflator, std::vector<Real> radiusBounds,
              std::vector<Real> translationBounds)
        : m_inflator(inflator), m_radiusBounds(radiusBounds),
          m_transBounds(translationBounds) { }

    template<class _Vector>
    void getParameterBounds(_Vector &lowerBounds, _Vector &upperBounds) {
        for (size_t p = 0; p < m_inflator.numParameters(); ++p) {
            switch (m_inflator.parameterType(p)) {
                case ParameterType::Thickness:
                    lowerBounds(p) = m_radiusBounds.at(0);
                    upperBounds(p) = m_radiusBounds.at(1);
                    break;
                case ParameterType::Offset:
                    lowerBounds(p) = m_transBounds.at(0);
                    upperBounds(p) = m_transBounds.at(1);
                    break;
                default: assert(false);
            }
        }
    }

    // Forward declaration so friend-ing can happen
    class IterationCallback;

    struct TensorFitCost : public ceres::CostFunction {
        typedef ceres::CostFunction Base;
        TensorFitCost(Inflator<_N> &inflator, const _ETensor &targetS)
            : m_inflator(inflator), m_targetS(targetS) {
            Base::set_num_residuals((_N == 2) ? 6 : 21);
            // We put all the pattern parameters in a single parameter block.
            Base::mutable_parameter_block_sizes()->assign(1,
                    m_inflator.numParameters());
        }

        virtual bool Evaluate(double const * const *parameters,
                double *residuals, double **jacobians) const {
            size_t nParams = parameter_block_sizes()[0];
            // Ceres seems to like re-evaluating at the same point, so we detect
            // this to avoid re-solving.
            if (!m_iterate || m_iterate->paramsDiffer(nParams, parameters[0])) {
                m_iterate = std::make_shared<Iterate>(m_inflator,
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
        Inflator<_N> &m_inflator;
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
                        nParams, &m_params[0], m_evaluator.m_targetS);
            }
            curr->writeDescription(std::cout);
            std::cout << std::endl;

            if (m_outPath != "")
                curr->writeMeshAndFields(m_outPath + "_" + std::to_string(m_iter));

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
        TensorFitCost *fitCost = new TensorFitCost(m_inflator, targetS);
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
            Iterate iterate(m_inflator, params.size(), params.data(), targetS);
            SField gradP = iterate.gradp_JS();
            iterate.writeDescription(std::cout);
            std::cout << std::endl;

            params -= gradP * alpha;
            // // TODO: remove this hack.
            // alpha *= .9;

            // Apply bound constraints
            params.minRelax(upperBounds);
            params.maxRelax(lowerBounds);

            if (outName != "")
                iterate.writeMeshAndFields(outName + "_" + std::to_string(i));
        }
    }

    // Evaluates the objective by inflating the wire mesh and homogenizing.
    // The iterate stored internally also knows how to evaluate the gradient
    // efficiently, so our GradientEvaluator below just accesses it.
    struct DLibObjectiveEvaluator {
        DLibObjectiveEvaluator(Inflator<_N> &inflator,
                const _ETensor &targetS)
            : m_iterate(NULL), m_inflator(inflator), m_targetS(targetS) {
            nParams = m_inflator.numParameters();
        };

        double operator()(const dlib_vector &x) const {
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);

            m_iterate = std::make_shared<Iterate>(m_inflator,
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
        Inflator<_N> &m_inflator;
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
                m_obj.currentIterate().writeMeshAndFields(m_outPath + "_" + std::to_string(m_iter));
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
        DLibObjectiveEvaluator obj(m_inflator, targetS);
        DLibGradientEvaluator grad(obj);

        size_t nParams = m_inflator.numParameters();
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
    Inflator<_N> &m_inflator;
    std::vector<Real> m_patternParams, m_radiusBounds, m_transBounds;
};

}

#endif /* end of include guard: PATTERNOPTIMIZATION_HH */
