// WARNING: THIS IS INCLUDED INSIDE THE "Optimizer" CLASS DEFINITION!


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

            // JP: unused...
			// Real JS  = curr->evaluateJS();
			// Real RT  = curr->evaluateRT(m_initialParams, m_weight);
			
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

	// MHS on Oct 2 2015:
	// this was written for use in the regularized version of the bfgs optimizer
	// TODO: see if you can update the original ReportingStopStartegy so that there
	// is no need for a new one !!!
	class ReportingStopStrategy2 : public dlib::objective_delta_stop_strategy {
        typedef dlib::objective_delta_stop_strategy Base;
    public:
        ReportingStopStrategy2(double min_delta, unsigned long max_iter,
                              DLibObjectiveEvaluator2 &obj, SField &initialParams, const double regularizationWeight, const std::string &outPath1, const std::string outPath2)
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
