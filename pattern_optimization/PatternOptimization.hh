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

#include <cassert>
#include <memory>
#include <map>

#include "Inflator.hh"
#include "PatternOptimizationIterate.hh"
#include <EdgeFields.hh>
#include <MSHFieldWriter.hh>

#include "PatternOptimizationConfig.hh"

// DLIB seems to interact badly with Boost if included early on; include it last.
#include <ceres/ceres.h>
#include <glog/logging.h>
#include <levmar.h>
#include <dlib/optimization.h>

namespace PatternOptimization {

}

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
    Optimizer(ConstrainedInflator<_N> &inflator, const std::vector<Real> &radiusBounds,
              const std::vector<Real> &translationBounds,
              const std::vector<Real> &blendingBounds,
              const std::map<size_t, Real> &varLowerBounds,
              const std::map<size_t, Real> &varUpperBounds)
        : m_inflator(inflator), m_radiusBounds(radiusBounds),
          m_transBounds(translationBounds), m_blendBounds(blendingBounds),
          m_varLowerBounds(varLowerBounds),
          m_varUpperBounds(varUpperBounds) { }

    template<class _Vector>
    void getParameterBounds(_Vector &lowerBounds, _Vector &upperBounds) {
        for (size_t p = 0; p < m_inflator.numParameters(); ++p) {
            // Explicitly specified bounds overried default type-based bounds.
            if (m_varLowerBounds.count(p)) {
                lowerBounds(p) = m_varLowerBounds.at(p);
                upperBounds(p) = m_varUpperBounds.at(p);
                continue;
            }
            switch (m_inflator.parameterType(p)) {
                case ParameterType::Thickness:
                    lowerBounds(p) = m_radiusBounds.at(0);
                    upperBounds(p) = m_radiusBounds.at(1);
                    break;
                case ParameterType::Offset:
                    lowerBounds(p) = m_transBounds.at(0);
                    upperBounds(p) = m_transBounds.at(1);
                    break;
                case ParameterType::Blending:
                    lowerBounds(p) = m_blendBounds.at(0);
                    upperBounds(p) = m_blendBounds.at(1);
                    break;
                default: assert(false);
            }
        }
    }

    // Enables creation and sharing of a iterates amongst several objects (e.g.
    // CeresCostWrapper(s) and IterationCallback)
    template<class _ItFactory>
    struct SharedIterateManager {
        SharedIterateManager(ConstrainedInflator<N> &inflator, const _ETensor &targetS)
            : m_inflator(inflator), m_targetS(targetS) { }

        Iterate &get(size_t nParams, double *params) {
            m_currIterate = _ItFactory::getIterate(std::move(m_currIterate),
                    m_inflator, nParams, parameters[0], m_targetS);
            return *m_currIterate;
        }

              Iterate &get()       { assert(m_currIterate); return *m_currIterate; }
        const Iterate &get() const { assert(m_currIterate); return *m_currIterate; }

        size_t numParameters() const { return m_inflator.numParameters(); }
        
    private:
        _ETensor m_targetS;
        ConstrainedInflator<_N> &m_inflator;
        std::unique_ptr<Iterate> m_currIterate;
    };

    // Wraps (nonlinear) least squares-supporting objective terms
    // (those with residuals and jacobians)
    template<class _ItManager>
    struct CeresCostWrapper : public ceres::CostFunction {
        typedef ceres::CostFunction Base;
        CeresCostWrapper(const std::string &quantity, _ItManager &itManager)
            : m_quantity(quantity), m_itManager(itManager)
        {
            const size_t nParams = m_itManager.numParameters();
            Base::set_num_residuals(Iterate::objectiveTermInfo(quantity, nParams).numResiduals());
            // We put all the pattern parameters in a single parameter block.
            Base::mutable_parameter_block_sizes()->assign(1, nParams);
        }

        virtual bool Evaluate(double const * const *parameters,
                double *residuals, double **jacobians) const {
            const size_t nParams = parameter_block_sizes()[0];
            auto &it = m_itManager.get(nParams, parameters[0]);

            const ObjectiveTerm &term = it.objectiveTerm(m_quantity);
            const size_t nResiduals = term.numResiduals();

            Real sqrt_weight = sqrt(term.normalizedWeight());
            for (size_t r = 0; r < nResiduals; ++r)
                residuals[r] = sqrt_weight * term.residual(r);

            if (jacobians == NULL) return true;

            for (size_t r = 0; r < nResiduals; ++r) {
                for (size_t p = 0; p < nParams; ++p)
                    jacobians[0][r * nParams + p] = sqrt_weight * term.jacobian(r, p);
            }

            return true;
        }

        const Iterate &currentIterate() const { return m_itManager.get(); }

        virtual ~CeresCostWrapper() { }

    private:
        const std::string m_quantity;
        mutable _ItManager &m_itManager;
    };

    template<class _ItManager>
    class IterationCallback : public ceres::IterationCallback {
    public:
        IterationCallback(_ItManager &itManager, SField &params,
                          const std::string &outPath)
            : m_itManager(itManager), m_params(params), m_outPath(outPath),
            m_iter(0) {}
        ceres::CallbackReturnType operator()(const ceres::IterationSummary &sum)
        {
            Iterate &curr = m_itManager.get(m_params.size(), &m_params[0]);
            curr->writeDescription(std::cout);
            std::cout << std::endl;

            if (m_outPath != "")
                curr->writeMeshAndFields(m_outPath + "_" + std::to_string(m_iter));

            ++m_iter;
            return ceres::SOLVER_CONTINUE;
        }

        virtual ~IterationCallback() { }
    private:
        _ItManager &m_itManager;
        SField &m_params;
        std::string m_outPath;
        size_t m_iter;
    };

    void optimize_lm(SField &params, const _ETensor &targetS,
                  const string &outPath) {
        using IM = SharedIterateManager<IterateFactory<>>;
        IM imanager(m_inflator, targetS);
        CeresCostWrapper<IM> *js = new CeresCostWrapper<IM>("JS", imanager);
        ceres::Problem problem;
        problem.AddResidualBlock(js, NULL, params.data());

        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        for (size_t p = 0; p < params.domainSize(); ++p) {
            problem.SetParameterLowerBound(params.data(), p, lowerBounds[p]);
            problem.SetParameterUpperBound(params.data(), p, upperBounds[p]);
        }

        ceres::Solver::Options options;
        options.update_state_every_iteration = true;
        IterationCallback<IM> cb(imanager, params, outPath);
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

    // Proximity-regularized optimization
	void optimize_lm_regularized(SField &params,
    							 SField &initialParams,
    							 const double regularizationWeight,
    							 const _ETensor &targetS,
    							 const string outPath, 
    							 Real & initialCost,
    							 Real & finalCost,
    							 _ETensor & stiffness) {
        using IM = SharedIterateManager<IterateFactory<IFConfigProximityRegularization>>;
        IM imanager(m_inflator, targetS);
        CeresCostWrapper<IM> *js  = new CeresCostWrapper<IM>("JS", imanager);
        CeresCostWrapper<IM> *reg = new CeresCostWrapper<IM>("ProximityRegularization", imanager);
        ceres::Problem problem;
        problem.AddResidualBlock(js, NULL, params.data());
        problem.AddResidualBlock(new Regularization(initialParams, regularizationWeight, params.domainSize()), NULL, params.data());

        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        for (size_t p = 0; p < params.domainSize(); ++p) {
            problem.SetParameterLowerBound(params.data(), p, lowerBounds[p]);
            problem.SetParameterUpperBound(params.data(), p, upperBounds[p]);
        }

        ceres::Solver::Options options;
        options.update_state_every_iteration = true;
        IterationCallback<IM> cb(imanager, params, outPath);
        options.callbacks.push_back(&cb);
        // options.minimizer_type = ceres::LINE_SEARCH;
        // options.line_search_direction_type = ceres::BFGS;
        // options.trust_region_strategy_type = ceres::DOGLEG;
        // options.dogleg_type = ceres::SUBSPACE_DOGLEG;
        // options.use_nonmonotonic_steps = true;
        // options.minimizer_progress_to_stdout = true;
        // options.initial_trust_region_radius = 1e4; // ceres's default
        options.max_num_iterations = 400;
        options.function_tolerance = 1.0e-16;
    	options.gradient_tolerance = 1.0e-32;
    	options.parameter_tolerance = 1.0e-32;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        /* std::cout << summary.BriefReport() << "\n"; */
        /* std::cout << summary.FullReport() << "\n"; */
        initialCost = summary.initial_cost;
        finalCost   = summary.final_cost;

		std::unique_ptr<Iterate> final_iter;
        Iterate &final_iter = imanager.get(params.size(), &params[0]);

		_ETensor final_C = final_iter.elasticityTensor();
		stiffness = final_C;

		std::cout << "the anisotopy of final C [computed in the optimizer before applying the transformations] is  " << final_C.anisotropy() << std::endl;
		std::cout << "testing final_iter to output JS " << final_iter.objectiveTerm("JS").evaluate() << std::endl;
		std::cout << "testing final_iter to output RT " << final_iter.objectiveTerm("ProximityRegularization").evaluate() << std::endl;

		std::cout << "this is the target C after the optimizer exits:" << std::endl;
		std::cout << "-----------------------------------------------" << std::endl;
		final_C.printOrthotropic(std::cout);
		std::cout << "-----------------------------------------------" << std::endl;

    }

    // MHS on AUG 25, 2015:
    // DOGLEG gets similar patterns (excluding rotations) for deformed cells with
    // Jacobian=[a 0; 0 b] and Jacobian=[b 0; 0 a]
    // TODO: more tests are needed to see if this is better than lm or not!
    void optimize_dogleg(SField &params, const _ETensor &targetS,
                  const string &outPath) {
        using IM = SharedIterateManager<IterateFactory<>>;
        IM imanager(m_inflator, targetS);
        CeresCostWrapper<IM> *js  = new CeresCostWrapper<IM>("JS", imanager);
        ceres::Problem problem;
        problem.AddResidualBlock(js, NULL, params.data());

        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        for (size_t p = 0; p < params.domainSize(); ++p) {
            problem.SetParameterLowerBound(params.data(), p, lowerBounds[p]);
            problem.SetParameterUpperBound(params.data(), p, upperBounds[p]);
        }

        ceres::Solver::Options options;
        options.update_state_every_iteration = true;
        IterationCallback<IM> cb(imanager, params, outPath);
        options.callbacks.push_back(&cb);
        // options.minimizer_type = ceres::LINE_SEARCH;
        // options.line_search_direction_type = ceres::BFGS;
        options.trust_region_strategy_type = ceres::DOGLEG;  // MHS: this seems to be needed in order
        options.dogleg_type = ceres::SUBSPACE_DOGLEG;        // to get similar (rotated) pattern
        options.use_nonmonotonic_steps = true;               // for F [a 0; 0 b] and [b 0; 0 a]
        // options.minimizer_progress_to_stdout = true;
        // options.initial_trust_region_radius = 0.01;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << "\n";
    }

    // Enforce bound constraints with a penalty block instead of
    // SetParameter{Upper,Lower}Bound
    struct BoundPenalty : public ceres::CostFunction {
        typedef ceres::CostFunction Base;
        BoundPenalty(size_t p, double lower, double upper, size_t numParams)
            : m_param(p), m_lower(lower), m_upper(upper) {
            assert(m_param < numParams);
            Base::set_num_residuals(1);
            // There is only a single block of pattern parameters
            Base::mutable_parameter_block_sizes()->assign(1, numParams);
        }

        // r = w * max(c^2 - 1, 0)
        // where c = (2p - u - l) / (u - l)
        virtual bool Evaluate(double const * const *parameters,
                double *residuals, double **jacobians) const {
            double p = parameters[0][m_param];
            size_t nParams = parameter_block_sizes()[0];

            double interval = m_upper - m_lower;
            double c = (2 * p - m_upper - m_lower) / interval;
            residuals[0] = m_weight * std::max(c * c - 1.0, 0.0);
            if (jacobians == NULL) return true;
            for (size_t i = 0; i < nParams; ++i) {
                jacobians[0][i] = 0.0;
            }
            if ((p < m_lower) || (p > m_upper))
                jacobians[0][m_param] = m_weight * (8 * p - 4 * (m_lower + m_upper)) / (interval * interval);
            return true;
        }

        virtual ~BoundPenalty() { }

    private:
        size_t m_param;
        double m_lower, m_upper;
        double m_weight = 1e10;
    };

    void optimize_lm_bound_penalty(SField &params, const _ETensor &targetS,
                  const string &outPath) {
        using IM = SharedIterateManager<IterateFactory<>>;
        IM imanager(m_inflator, targetS);
        CeresCostWrapper<IM> *js  = new CeresCostWrapper<IM>("JS", imanager);
        ceres::Problem problem;
        problem.AddResidualBlock(js, NULL, params.data());

        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        for (size_t p = 0; p < params.domainSize(); ++p) {
            problem.AddResidualBlock(
                    new BoundPenalty(p, lowerBounds[p], upperBounds[p], nParams),
                    NULL, params.data());
        }

        ceres::Solver::Options options;
        options.update_state_every_iteration = true;
        IterationCallback<IM> cb(imanager, params, outPath);
        options.callbacks.push_back(&cb);

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << "\n";
    }

    template<class _ItManager>
    struct LevmarEvaluator {
        LevmarEvaluator(_ItManager &imanager) : m_imanager(imanager) { }

        void residual(double *params, double *residual, int numParams, int numResiduals) {
            auto &it = m_imanager.get(numParams, params);
            size_t offset = 0;
            for (const auto &term : it.objectiveTerms()) {
                auto &nllsTerm = dynamic_cast<NLLSObjectiveTerm&>(*term.second);
                const size_t nr = Iterate::objectiveTermInfo(term.first, numParams).numResiduals();
                Real sqrt_weight = sqrt(term.normalizedWeight());

                for (size_t r = 0; r < nr; ++r)
                    residuals[offset + r] = sqrt_weight * nllsTerm.residual(r);
                offset += nr;
            }

            assert(offset == (size_t) numResiduals);
        }
    
        // Row major Jacobian
        void jacobian(double *params, double *jacobian, int numParams, int numResiduals) {
            auto &it = m_imanager.get(numParams, params);
            size_t offset = 0;
            for (const auto &term : it.objectiveTerms()) {
                auto &nllsTerm = dynamic_cast<NLLSObjectiveTerm&>(*term.second);
                const size_t nr = Iterate::objectiveTermInfo(term.first, numParams).numResiduals();
                Real sqrt_weight = sqrt(term.normalizedWeight());

                for (size_t r = 0; r < nr; ++r) {
                    for (size_t p = 0; p < size_t(numParams); ++p)
                        jacobian[(offset + r) * numParams + p] = sqrt_weight * nllsTerm.jacobian(r, p);
                }
                offset += nr;
            }

            assert(offset == (size_t) numResiduals);
        }

        // Residual and Jacobian evaluation callbacks for levmar
        // The Jacobian is stored row-major.
        static void residual(double *params, double *residual, int numParams, int numResiduals, void *instance) { static_cast<LevmarEvaluator *>(instance)->residual(params, residual, numParams, numResiduals); }
        static void jacobian(double *params, double *jacobian, int numParams, int numResiduals, void *instance) { static_cast<LevmarEvaluator *>(instance)->jacobian(params, jacobian, numParams, numResiduals); }

    private:
        _ItManager &m_imanager;
    };

    void optimize_levmar(SField &params, const _ETensor &targetS, size_t niters, const string &outName) {
        size_t numResiduals = ObjectiveTermInfo::get<_N>(nParams).numResiduals();
        size_t numParams = params.domainSize();
        SField lowerBounds(numParams), upperBounds(numParams);
        getParameterBounds(lowerBounds, upperBounds);

        using IM = SharedIterateManager<IterateFactory<>>;
        IM imanager(m_inflator, targetS);
        LevmarEvaluator<IM> eval(imanager);

        double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
        opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
        opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 
        // TODO: add iteration completed callback.
        dlevmar_bc_der(LevmarEvaluator<IM>::residual, LevmarEvaluator<IM>::jacobian,
                       params.data(),
                       NULL, // Since our "measurements" are actually residuals, target x is 0
                       numParams, numResiduals,
                       lowerBounds.data(), upperBounds.data(),
                       NULL, // dscl--no diagonal scaling
                       niters, opts, info, NULL, NULL, static_cast<void *>(&eval));
        for (size_t i = 0; i < LM_INFO_SZ; ++i)
            cout << info[i] << "\t";
    }

    void optimize_gd(SField &params, const _ETensor &targetS,
            size_t niters, double alpha, const string &outName) {
        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        
        using IM = SharedIterateManager<IterateFactory<>>;
        IM imanager(m_inflator, targetS);

        // Gradient descent version
        // TODO: switch to steepestDescentParam when implemented
        for (size_t i = 0; i < niters; ++i) {
            auto &it = imanager.get(params.size(), params.data());
            SField gradp = it.gradp(m_inflator.shapeVelocities());
            it.writeDescription(std::cout);
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
    template<class _ItManager>
    struct DLibObjectiveEvaluator {
        DLibObjectiveEvaluator(_ItManager &imanager) : m_imanager(imanager) {
            nParams = m_imanager.numParameters();
        };

        double operator()(const dlib_vector &x) const {
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);

            auto &it = m_imanager.get(nParams, &x_vec[0]);
            return it.evaluate();
        }

        const Iterate &currentIterate() const { return m_imanager.get(); }

        size_t nParams;
    private:
        // Iterate manager is mutable so that operator() can be const as dlib requires
        mutable _ItManager m_imanager;
    };
    
    // Extracts gradient from the iterate constructed by DLibObjectiveEvaluator.
    template<class _ItManager>
    struct DLibGradientEvaluator {
        DLibGradientEvaluator(const DLibObjectiveEvaluator<_ItManager> &obj) : m_obj(obj) { }

        dlib_vector operator()(const dlib_vector &x) const {
            size_t nParams = m_obj.nParams;
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);
            if (m_obj.currentIterate().paramsDiffer(nParams, &x_vec[0]))
                throw std::runtime_error("Objective must be evaluated first");
            SField gradp(m_obj.currentIterate().gradp());

            dlib_vector result(nParams);
            for (size_t p = 0; p < nParams; ++p)
                result(p) = gradp[p];
            return result;
        }
    private:
        const DLibObjectiveEvaluator<_ItManager> &m_obj;
    };

    // Hack to get notified at the end of each iteration: subclass the stop
    // strategy.
    template<class _ItManager>
    class ReportingStopStrategy : public dlib::objective_delta_stop_strategy {
        typedef dlib::objective_delta_stop_strategy Base;
    public:
        ReportingStopStrategy(double min_delta, unsigned long max_iter,
                              DLibObjectiveEvaluator<_ItManager> &obj, const std::string &outPath)
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
        const DLibObjectiveEvaluator<_ItManager> &m_obj;
        std::string m_outPath;
        size_t m_iter;
    };

    // If max_size = 0, plain bfgs is used
    // otherwise l-bfgs is used.
    void optimize_bfgs(SField &params, const _ETensor &targetS, size_t niters,
                       const string &outPath, size_t max_size = 0) {
        using IM = SharedIterateManager<IterateFactory<>>;
        IM imanager(m_inflator, targetS);

        DLibObjectiveEvaluator<IM> obj(imanager);
        DLibGradientEvaluator<IM> grad(obj);

        size_t nParams = m_inflator.numParameters();
        // convert initial parameter vector
        dlib_vector optParams(nParams);
        for (size_t p = 0; p < nParams; ++p)
            optParams(p) = params[p];

        dlib_vector lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);

        if (max_size == 0)
            dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
                    ReportingStopStrategy<IM>(1e-16, niters, obj, outPath),
                    obj, grad, optParams, lowerBounds, upperBounds);
        else
            dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(max_size),
                    ReportingStopStrategy<IM>(1e-16, niters, obj, outPath),
                    obj, grad, optParams, lowerBounds, upperBounds);

        // convert solution
        for (size_t p = 0; p < nParams; ++p)
            params[p] = optParams(p);
    }


private:
    ConstrainedInflator<_N> &m_inflator;
    std::vector<Real> m_patternParams, m_radiusBounds, m_transBounds, m_blendBounds;
    std::map<size_t, Real> m_varLowerBounds, m_varUpperBounds;

}; // class Optimizer

} // namespace PatternOptimization

#endif /* end of include guard: PATTERNOPTIMIZATION_HH */
