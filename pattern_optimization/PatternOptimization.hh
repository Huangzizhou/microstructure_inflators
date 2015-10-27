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
#include <levmar.h>

#include <cassert>
#include <memory>
#include <map>

#include "Inflator.hh"
#include "PatternOptimizationIterate.hh"
#include <EdgeFields.hh>
#include <MSHFieldWriter.hh>

#include "PatternOptimizationConfig.hh"

namespace PatternOptimization {

// Use previous iterate if evaluating the same point. Otherwise, attempt to
// inflate the new parameters. Try three times to inflate, and if unsuccessful,
// estimate the point with linear extrapolation.
template<class _Iterate, class _Inflator, class _ETensor>
std::shared_ptr<_Iterate>
getIterate(std::shared_ptr<_Iterate> oldIterate,
        _Inflator &inflator, size_t nParams, const double *params,
        const _ETensor &targetS) {
    if (!oldIterate || oldIterate->paramsDiffer(nParams, params)) {
        std::shared_ptr<_Iterate> newIterate;
        bool success;
        for (size_t i = 0; i < 3; ++i) {
            success = true;
            try {
                newIterate = std::make_shared<_Iterate>(inflator,
                                nParams, params, targetS);
            }
            catch (std::exception &e) {
                std::cerr << "INFLATOR FAILED: " << e.what() << endl;
                success = false;
            }
            if (success) break;
        }
        if (!success) {
            std::cerr << "3 INFLATION ATTEMPTS FAILED." << std::endl;
            if (!oldIterate) throw std::runtime_error("Inflation failure on first iterate");
            newIterate = oldIterate;
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
              const std::map<size_t, Real> &varLowerBounds,
              const std::map<size_t, Real> &varUpperBounds)
        : m_inflator(inflator), m_radiusBounds(radiusBounds),
          m_transBounds(translationBounds), m_varLowerBounds(varLowerBounds),
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
                default: assert(false);
            }
        }
    }

    // Forward declaration so friend-ing can happen
    class IterationCallback;
    class IterationCallback2;
    class IterationCallback3;

    struct TensorFitCost : public ceres::CostFunction {
        typedef ceres::CostFunction Base;
        TensorFitCost(ConstrainedInflator<_N> &inflator, const _ETensor &targetS)
            : m_inflator(inflator), m_targetS(targetS) {
            Base::set_num_residuals((_N == 2) ? 6 : 21);
            // We put all the pattern parameters in a single parameter block.
            Base::mutable_parameter_block_sizes()->assign(1,
                    m_inflator.numParameters());
        }

        virtual bool Evaluate(double const * const *parameters,
                double *residuals, double **jacobians) const {
            size_t nParams = parameter_block_sizes()[0];
            m_iterate = getIterate(m_iterate, m_inflator, nParams,
                                   parameters[0], m_targetS);

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
        ConstrainedInflator<_N> &m_inflator;
        _ETensor m_targetS;
        // Ceres requires Evaluate to be constant, so this caching pointer must
        // be made mutable.
        mutable std::shared_ptr<Iterate> m_iterate;

        friend class IterationCallback;
        friend class IterationCallback2;
        friend class IterationCallback3;
    };

    class IterationCallback : public ceres::IterationCallback {
    public:
        IterationCallback(TensorFitCost &evalulator, SField &params,
                const std::string &outPath)
            : m_evaluator(evalulator), m_params(params), m_outPath(outPath),
            m_iter(0) {}
        ceres::CallbackReturnType operator()(const ceres::IterationSummary &sum)
        {
            auto curr = getIterate(m_evaluator.m_iterate,
                            m_evaluator.m_inflator, m_params.size(), &m_params[0],
                            m_evaluator.m_targetS);
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

	// MHS on Oct 2 2015:
	// this was written for use with lm_regularized
	// TODO: see if you can avoid it by using the original "IterationCallback"
   	class IterationCallback2 : public ceres::IterationCallback {
    public:
        IterationCallback2(TensorFitCost &evalulator, SField &initialParams, SField &params, const double regularizationWeight,
                const std::string &outPath1, const std::string &outPath2)
            : m_evaluator(evalulator), m_initialParams(initialParams), m_params(params), m_weight(regularizationWeight), m_outPath1(outPath1), m_outPath2(outPath2),
            m_iter(0) {}
        ceres::CallbackReturnType operator()(const ceres::IterationSummary &sum)
        {
        	auto curr = getIterate(m_evaluator.m_iterate,
        			m_evaluator.m_inflator, m_params.size(), &m_params[0],
        			m_evaluator.m_targetS);

			curr->writeDescription(std::cout);
			std::cout << std::endl;

			if (m_outPath2 != "")
				curr->writeMeshAndFields(m_outPath2 + "/_" + std::to_string(m_iter));

			double JS  = curr->evaluateJS();
			double RT  = curr->evaluateRT(m_initialParams, m_weight);
			auto gradp_JS = curr->gradp_JS();
			auto gradp_RT = curr->gradp_RT(m_initialParams, m_weight);

			std::vector<double> F;
			std::vector<std::vector<double>> J;

			F.clear();
			J.clear();

			for (size_t i = 0; i < 3; ++i)
				for (size_t j = i; j < 3; ++j)
					F.push_back(curr->residual(i, j));

            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = i; j < 3; ++j) {
					std::vector<double> jacobianRow;
					jacobianRow.clear();
                    for (size_t p = 0; p < m_params.size(); ++p)
                        jacobianRow.push_back(curr->jacobian(i, j, p));
                    J.push_back(jacobianRow);
                }
            }

			std::cout << "----- Jacobian Matrix -----" << std::endl;
			for (auto jRow = J.begin(); jRow != J.end(); ++jRow){
				for (auto jElm = (*jRow).begin(); jElm != (*jRow).end(); ++jElm)
					std::cout << *jElm << "\t";
				std::cout << std::endl;
			}
			std::cout << "---------------------------" << std::endl;

			std::cout << "----- Residual Vector -----" << std::endl;
			for (auto fElm = F.begin(); fElm != F.end(); ++fElm)
				std::cout << *fElm << std::endl;
			std::cout << "---------------------------" << std::endl;


			std::cout << "-------" << std::endl;
			gradp_JS.print(std::cout, "", "", "", "\t");
			std::cout << std::endl << gradp_JS.norm() << std::endl << sum.gradient_norm << std::endl;

			/* int i; */
			/* i = getchar(); */

			std::ofstream ofs1, ofs2;

			ofs1.open(m_outPath1 + "/iterations_optimization.txt", std::ofstream::out | std::ofstream::app);
			ofs2.open(m_outPath1 + "/iterations_parameters.txt",   std::ofstream::out | std::ofstream::app);

			ofs1.close();
			ofs2.close();

			int flag1 = 0;
			int flag2 = 0;
			std::string line;
			std::ifstream ifs1, ifs2;
			ifs1.open(m_outPath1 + "/iterations_optimization.txt");
			ifs2.open(m_outPath1 + "/iterations_parameters.txt");
			if (!getline(ifs1, line))
				flag1 = 1;
			if (!getline(ifs2, line))
				flag2 = 1;
			ifs1.close();
			ifs1.close();


			ofs1.open(m_outPath1 + "/iterations_optimization.txt", std::ofstream::out | std::ofstream::app);
			ofs2.open(m_outPath1 + "/iterations_parameters.txt",   std::ofstream::out | std::ofstream::app);

			if (flag1 == 1)
				ofs1 << "JS" << "\t" << "RT" << "\t" << "grad\\_JS" << "\t" << "grad\\_RT" << "\t" <<
					    "Iter" << "\t" << "Cost" << "\t" << "grad_norm" << "\t" << "rel_decrease" << "\t" << "cost_change" << "\t" <<  "trust_radius" << std::endl;

			if (flag2 == 1)
				ofs2 << "p1" << "\t" << "p2" << "\t" << "p3" << "\t" << "p4" << "\t" << "p5" << "\t" << "p6" << "\t" << "p7" << "\t" << "p8" << "\t" << "p9" << std::endl;


			ofs1 << std::setprecision(16) << std::showpos;
			ofs2 << std::setprecision(16) << std::showpos;


			if (m_iter == 0)
			{
				ofs1 << ">DATASET<" << std::endl;
				ofs2 << ">DATASET<" << std::endl;
			}


			ofs1 << JS << "\t" << RT << "\t" << gradp_JS.norm() << "\t" << gradp_RT.norm() << "\t" <<
					sum.iteration << "\t" << sum.cost << "\t" << sum.gradient_norm << "\t" << sum.relative_decrease << "\t" << sum.cost_change << "\t" <<  sum.trust_region_radius << std::endl;
			ofs1.close();

			for (size_t i = 0; i < m_params.size(); ++i)
				ofs2 << m_params[i] << "\t";
			ofs2 << std::endl;

            ++m_iter;
            return ceres::SOLVER_CONTINUE;
        }

        virtual ~IterationCallback2() { }
    private:
        TensorFitCost &m_evaluator;
        SField &m_params;
        SField &m_initialParams;
        double m_weight;
        std::string m_outPath1, m_outPath2;
        size_t m_iter;
        size_t m_flag;
    };

	// MHS on Oct 26 2015:
	// this is a new callback for use with lm_regularized
	// TODO: see if you can avoid it by using the original "IterationCallback"
   	class IterationCallback3 : public ceres::IterationCallback {
    public:
        IterationCallback3(TensorFitCost &evalulator, SField &initialParams, SField &params, 
						   _ETensor &stiffness, 
						   const double regularizationWeight,
						   const string outPath)
            : m_evaluator(evalulator), m_initialParams(initialParams), m_params(params), m_weight(regularizationWeight), m_outPath(outPath), m_stiffness(&stiffness) {}

        ceres::CallbackReturnType operator()(const ceres::IterationSummary &sum)
        {
        	auto curr = getIterate(m_evaluator.m_iterate,
        			m_evaluator.m_inflator, m_params.size(), &m_params[0],
        			m_evaluator.m_targetS);

            if (m_outPath != "")
                curr->writeMeshAndFields(m_outPath + "_" + std::to_string(m_iter));

			curr->writeDescription(std::cout);
			std::cout << std::endl;

			Real JS  = curr->evaluateJS();
			Real RT  = curr->evaluateRT(m_initialParams, m_weight);
			
			_ETensor C = curr->elasticityTensor();
			
			*m_stiffness      = C;

            ++m_iter;

            return ceres::SOLVER_CONTINUE;
        }
        /* Real &getInitialCost()   {return m_initialCost;} */
        /* Real &getTotalCost()     {return m_totalCost;} */
        /* Real &getStifnessCost()  {return m_stiffnessCost;} */
        /* _ETensor &getStiffness() {return m_stiffness;} */

        virtual ~IterationCallback3() { }
    private:
        TensorFitCost &m_evaluator;
        SField &m_params;
        SField &m_initialParams;
        double m_weight;
        size_t m_iter;
        string m_outPath;
        _ETensor * m_stiffness;
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


    // MHS on AUG 31, 2015:
	// adding the following  regularization to the cost function
	// alpha * (p - p0).(p - p0) / 2
    struct Regularization : public ceres::CostFunction {
        typedef ceres::CostFunction Base;
        Regularization(const SField & initialParams, const double regularizationWeight, size_t numParams)
            : m_initialParams(initialParams), m_weight(regularizationWeight) {
            Base::set_num_residuals(numParams);
            // There is only a single block of pattern parameters
            Base::mutable_parameter_block_sizes()->assign(1, numParams);
        }

        virtual bool Evaluate(double const * const *parameters,
                double *residuals, double **jacobians) const {

            size_t nParams = parameter_block_sizes()[0];

			for (size_t i = 0; i < nParams; ++i)
				residuals[i] = std::sqrt(m_weight) * (parameters[0][i] - m_initialParams[i]);

            if (jacobians == NULL) return true;
            for (size_t i = 0; i < nParams; ++i)
            	for (size_t j = 0; j < nParams; ++j)
                	jacobians[0][i * nParams + j] = (i == j) ? std::sqrt(m_weight) : 0.0;

            return true;
        }

        virtual ~Regularization() { }

    private:
        SField m_initialParams;
        double m_weight;
    };

	// MHS on Oct 2:
	// the regularized lm optimizer
	void optimize_lm_regularized(SField &params,
    							 SField &initialParams,
    							 const double regularizationWeight,
    							 const _ETensor &targetS,
    							 const string outPath, 
    							 Real & initialCost,
    							 Real & finalCost,
    							 _ETensor & stiffness) {
        TensorFitCost *fitCost = new TensorFitCost(m_inflator, targetS);
        ceres::Problem problem;
        problem.AddResidualBlock(fitCost, NULL, params.data());
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
        IterationCallback3 cb(*fitCost, initialParams, params, stiffness, regularizationWeight, outPath);
        options.callbacks.push_back(&cb);
        // options.minimizer_type = ceres::LINE_SEARCH;
        // options.line_search_direction_type = ceres::BFGS;
        // options.trust_region_strategy_type = ceres::DOGLEG;
        // options.dogleg_type = ceres::SUBSPACE_DOGLEG;
        // options.use_nonmonotonic_steps = true;
        // options.minimizer_progress_to_stdout = true;
        // options.initial_trust_region_radius = 1e4; // ceres's default
        options.max_num_iterations = 100;
        /* options.function_tolerance = 1.0e-16; */
    	/* options.gradient_tolerance = 1.0e-32; */
    	/* options.parameter_tolerance = 1.0e-32; */
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << "\n";
        std::cout << summary.FullReport() << "\n";
        initialCost = summary.initial_cost;
        finalCost   = summary.final_cost;
    }

    // MHS on AUG 25, 2015:
    // DOGLEG gets similar patterns (excluding rotations) for deformed cells with
    // Jacobian=[a 0; 0 b] and Jacobian=[b 0; 0 a]
    // TODO: more tests are needed to see if this is better than lm or not!
    void optimize_dogleg(SField &params, const _ETensor &targetS,
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
        TensorFitCost *fitCost = new TensorFitCost(m_inflator, targetS);
        ceres::Problem problem;
        problem.AddResidualBlock(fitCost, NULL, params.data());

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
        IterationCallback cb(*fitCost, params, outPath);
        options.callbacks.push_back(&cb);

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << "\n";
    }

    struct LevmarEvaluator {
        LevmarEvaluator(ConstrainedInflator<_N> &inflator, const _ETensor &targetS)
            : m_inflator(inflator), m_targetS(targetS) { }

        void residual(double *params, double *residual, int numParams, int numResiduals) {
            m_iterate = getIterate(m_iterate, m_inflator, numParams, params, m_targetS);
            size_t r = 0;
            for (size_t i = 0; i < flatLen(_N); ++i)
                for (size_t j = i; j < flatLen(_N); ++j)
                    residual[r++] = m_iterate->residual(i, j);
            assert(r == (size_t) numResiduals);
        }
    
        // Row major Jacobian
        void jacobian(double *params, double *jacobian, int numParams, int numResiduals) {
            m_iterate = getIterate(m_iterate, m_inflator, numParams, params, m_targetS);
            size_t r = 0;
            for (size_t i = 0; i < flatLen(_N); ++i) {
                for (size_t j = i; j < flatLen(_N); ++j) {
                    for (size_t p = 0; p < numParams; ++p)
                        jacobian[r * numParams + p] = m_iterate->jacobian(i, j, p);
                    ++r;
                }
            }
            assert(r == (size_t) numResiduals);
        }

        // Residual and Jacobian evaluation callbacks for levmar
        // The Jacobian is stored row-major.
        static void residual(double *params, double *residual, int numParams, int numResiduals, void *instance) { static_cast<LevmarEvaluator *>(instance)->residual(params, residual, numParams, numResiduals); }
        static void jacobian(double *params, double *jacobian, int numParams, int numResiduals, void *instance) { static_cast<LevmarEvaluator *>(instance)->jacobian(params, jacobian, numParams, numResiduals); }

        ConstrainedInflator<_N> &m_inflator;
        _ETensor m_targetS;
        std::shared_ptr<Iterate> m_iterate;
    };

    void optimize_levmar(SField &params, const _ETensor &targetS, size_t niters, const string &outName) {
        size_t numResiduals = (_N == 2) ? 6 : 21;
        size_t numParams = params.domainSize();
        SField lowerBounds(numParams), upperBounds(numParams);
        getParameterBounds(lowerBounds, upperBounds);

        LevmarEvaluator eval(m_inflator, targetS);

        double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
        opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
        opts[4]= LM_DIFF_DELTA; // relevant only if the Jacobian is approximated using finite differences; specifies forward differencing 
        // TODO: add iteration completed callback.
        dlevmar_bc_der(LevmarEvaluator::residual, LevmarEvaluator::jacobian,
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
        DLibObjectiveEvaluator(ConstrainedInflator<_N> &inflator,
                const _ETensor &targetS)
            : m_iterate(NULL), m_inflator(inflator), m_targetS(targetS) {
            nParams = m_inflator.numParameters();
        };

        double operator()(const dlib_vector &x) const {
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);

            m_iterate = getIterate(m_iterate, m_inflator, nParams, &x_vec[0],
                                   m_targetS);
            return m_iterate->evaluateJS();
        }

        const Iterate &currentIterate() const {
            assert(m_iterate); return *m_iterate;
        };

        size_t nParams;
    private:
        // Iterate is mutable so that operator() can be const as dlib requires
        mutable std::shared_ptr<Iterate> m_iterate;
        ConstrainedInflator<_N> &m_inflator;
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

  	// MHS on Oct 2 2015:
	// the following two structs were written for use with the regularized version of the bfgs optimizer
	// TODO: see if you can update the original "DLibObjectiveEvaluator" and "DLibGradientEvaluator" i
	// so that there is no need for the new ones !!!
  	struct DLibObjectiveEvaluator2 {
        DLibObjectiveEvaluator2(ConstrainedInflator<_N> &inflator,
                const _ETensor &targetS, SField &initialParams, const double regularizationWeight)
            : m_iterate(NULL), m_inflator(inflator), m_targetS(targetS), m_initialParams(initialParams), m_weight(regularizationWeight) {
            nParams = m_inflator.numParameters();
        };

        double operator()(const dlib_vector &x) const {
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);

            m_iterate = getIterate(m_iterate, m_inflator, nParams, &x_vec[0],
                                   m_targetS);
            return m_iterate->evaluateJS() + m_iterate->evaluateRT(m_initialParams, m_weight);
        }

        const Iterate &currentIterate() const {
            assert(m_iterate); return *m_iterate;
        };

        size_t nParams;
    private:
        // Iterate is mutable so that operator() can be const as dlib requires
        mutable std::shared_ptr<Iterate> m_iterate;
        ConstrainedInflator<_N> &m_inflator;
        _ETensor m_targetS;
        SField &m_initialParams;
        double m_weight;
    };

    struct DLibGradientEvaluator2 {
        DLibGradientEvaluator2(const DLibObjectiveEvaluator2 &obj, SField &initialParams, const double regularizationWeight)
        	: m_obj(obj), m_initialParams(initialParams), m_weight(regularizationWeight) { }

        dlib_vector operator()(const dlib_vector &x) const {
            size_t nParams = m_obj.nParams;
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);
            if (m_obj.currentIterate().paramsDiffer(nParams, &x_vec[0]))
                throw std::runtime_error("Objective must be evaluated first");
            SField gradp_JS(m_obj.currentIterate().gradp_JS());
            SField gradp_RT(m_obj.currentIterate().gradp_RT(m_initialParams, m_weight));

            dlib_vector result(nParams);
            for (size_t p = 0; p < nParams; ++p)
                result(p) = gradp_JS[p] + gradp_RT[p];
            return result;
        }
    private:
        const DLibObjectiveEvaluator2 &m_obj;
        double m_weight;
        SField &m_initialParams;
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


	// MHS on Oct 2 2015:
	// this was written for use in the regularized version of the bfgs optimizer
	// TODO: see if you can update the original ReportingStopStartegy so that there
	// is no need for a new one !!!
	class ReportingStopStrategy2 : public dlib::objective_delta_stop_strategy {
        typedef dlib::objective_delta_stop_strategy Base;
    public:
        ReportingStopStrategy2(double min_delta, unsigned long max_iter,
                              DLibObjectiveEvaluator2 &obj, SField initialParams, const double regularizationWeight, const std::string &outPath1, const std::string outPath2)
            : Base(min_delta, max_iter), m_obj(obj), m_outPath1(outPath1), m_outPath2(outPath2), m_initialParams(initialParams), m_weight(regularizationWeight), m_iter(0) { }

        template <typename T>
        bool should_continue_search(const T& x, const double funct_value,
            const T& funct_derivative) {
            m_obj.currentIterate().writeDescription(std::cout);
            cout << endl;

			auto curr = m_obj.currentIterate();


            if (m_outPath2 != "")
                curr.writeMeshAndFields(m_outPath2 + "/_" + std::to_string(m_iter));

			double JS  = curr.evaluateJS();
			double RT  = curr.evaluateRT(m_initialParams, m_weight);
			auto gradp_JS = curr.gradp_JS();
			auto gradp_RT = curr.gradp_RT(m_initialParams, m_weight);
			SField currentParams = curr.getParams();

			std::ofstream ofs1, ofs2;

			ofs1.open(m_outPath1 + "/iterations_optimization.txt", std::ofstream::out | std::ofstream::app);
			ofs2.open(m_outPath1 + "/iterations_parameters.txt",   std::ofstream::out | std::ofstream::app);

			ofs1.close();
			ofs2.close();

			int flag1 = 0;
			int flag2 = 0;
			std::string line;
			std::ifstream ifs1, ifs2;
			ifs1.open(m_outPath1 + "/iterations_optimization.txt");
			ifs2.open(m_outPath1 + "/iterations_parameters.txt");
			if (!getline(ifs1, line))
				flag1 = 1;
			if (!getline(ifs2, line))
				flag2 = 1;
			ifs1.close();
			ifs1.close();


			ofs1.open(m_outPath1 + "/iterations_optimization.txt", std::ofstream::out | std::ofstream::app);
			ofs2.open(m_outPath1 + "/iterations_parameters.txt",   std::ofstream::out | std::ofstream::app);

			if (flag1 == 1)
				ofs1 << "grad\\_JS" << "\t" << "grad\\_RT" << "\t" << "JS" << "\t" << "RT" << std::endl;

			if (flag2 == 1)
				ofs2 << "p1" << "\t" << "p2" << "\t" << "p3" << "\t" << "p4" << "\t" << "p5" << "\t" << "p6" << "\t" << "p7" << "\t" << "p8" << "\t" << "p9" << std::endl;


			ofs1 << std::setprecision(16) << std::showpos;
			ofs2 << std::setprecision(16) << std::showpos;

			if (m_iter == 0)
			{
				ofs1 << ">DATASET<" << std::endl;
				ofs2 << ">DATASET<" << std::endl;
			}


			ofs1 << JS << "\t" << RT << "\t" << gradp_JS.norm() << "\t" << gradp_RT.norm() << std::endl;
			ofs1.close();

			for (size_t i = 0; i < currentParams.size(); ++i)
				ofs2 << currentParams[i] << "\t";
			ofs2 << std::endl;

            ++m_iter;
            return Base::should_continue_search(x, funct_value, funct_derivative);
        }

    private:
        const DLibObjectiveEvaluator2 &m_obj;
        std::string m_outPath1, m_outPath2;
        size_t m_iter;
        double m_weight;
        SField &m_initialParams;
    };


	// MHS on Oct 2 2015:
	// TODO: this could be implemented better
	// TODO: more tests need to be run on this regularized version
	void optimize_bfgs_regularized(SField &params,
								   SField &initialParams,
   		                           const double regularizationWeight,
   		                           const _ETensor &targetS,
                       			   const string &outPath1,
                       			   const string &outPath2,
   		                           size_t niters,
                       			   size_t max_size = 0) {

        DLibObjectiveEvaluator2 obj(m_inflator, targetS, initialParams, regularizationWeight);
        DLibGradientEvaluator2 grad(obj, initialParams, regularizationWeight);

        size_t nParams = m_inflator.numParameters();
        // convert initial parameter vector
        dlib_vector optParams(nParams);
        for (size_t p = 0; p < nParams; ++p)
            optParams(p) = params[p];

        dlib_vector lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);

        if (max_size == 0)
            dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
                    ReportingStopStrategy2(1e-16, niters, obj, initialParams, regularizationWeight, outPath1, outPath2),
                    obj, grad, optParams, lowerBounds, upperBounds);
        else
            dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(max_size),
                    ReportingStopStrategy2(1e-16, niters, obj, initialParams, regularizationWeight, outPath1, outPath2),
                    obj, grad, optParams, lowerBounds, upperBounds);

        // convert solution
        for (size_t p = 0; p < nParams; ++p)
            params[p] = optParams(p);
    }
private:
    ConstrainedInflator<_N> &m_inflator;
    std::vector<Real> m_patternParams, m_radiusBounds, m_transBounds;
    std::map<size_t, Real> m_varLowerBounds, m_varUpperBounds;
};

}

#endif /* end of include guard: PATTERNOPTIMIZATION_HH */
