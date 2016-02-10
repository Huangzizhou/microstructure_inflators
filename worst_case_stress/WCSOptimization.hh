////////////////////////////////////////////////////////////////////////////////
// WCSOptimization.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Pattern optimization code for worst-case stress minimization (and
//      homogenized tensor fitting)
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/16/2015 13:12:44
////////////////////////////////////////////////////////////////////////////////
#ifndef WCSOPTIMIZATION_HH
#define WCSOPTIMIZATION_HH

#include <dlib/optimization.h>

#include <Inflator.hh>
#include <MSHFieldWriter.hh>

#include <ElasticityTensor.hh>

#include "WCSObjective.hh"
#include <PatternOptimizationConfig.hh>
#include <PatternOptimizationIterate.hh>

namespace WCStressOptimization {

template<class _Sim, class _Inflator, template<class> class _Iterate>
class Optimizer
{
    typedef typename _Sim::VField VField;
    typedef ScalarField<Real> SField;
    typedef typename _Sim::ETensor _ETensor;
    typedef dlib::matrix<double,0,1> dlib_vector;
    static constexpr size_t _N = _Sim::N;
public:
    Optimizer(_Inflator &inflator, const std::vector<Real> &radiusBounds,
              const std::vector<Real> &translationBounds,
              const std::vector<Real> &blendingBounds,
              const std::map<size_t, Real> &varLowerBounds,
              const std::map<size_t, Real> &varUpperBounds)
        : m_inflator(inflator), m_radiusBounds(radiusBounds),
          m_transBounds(translationBounds), m_blendingBounds(blendingBounds),
          m_varLowerBounds(varLowerBounds), m_varUpperBounds(varUpperBounds) { }

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
                    lowerBounds(p) = m_blendingBounds.at(0);
                    upperBounds(p) = m_blendingBounds.at(1);
                    break;
                default: assert(false);
            }
        }
    }
    
    void optimize_gd(SField &params, WCStressOptimization::Objective<_N> &fullObjective,
            size_t niters, double stepSize, const string &outName) {
        size_t nParams = params.domainSize();
        SField lowerBounds(nParams), upperBounds(nParams);
        getParameterBounds(lowerBounds, upperBounds);
        // Gradient Descent Version
        for (size_t i = 0; i < niters; ++i) {
            _Iterate<_Sim> iterate(m_inflator, params.size(), params.data(), fullObjective);
            iterate.writeDescription(std::cout);
            std::cout << std::endl;

            SField gradp_JFull = iterate.gradp_JFull();

            params -= gradp_JFull * stepSize;

            // Apply bound constraints
            params.minRelax(upperBounds);
            params.maxRelax(lowerBounds);

            if (outName != "")
                iterate.writeMeshAndFields(outName + "_" + std::to_string(i) + ".msh");
        }
    }

    // Evaluates the objective by inflating the wire mesh and homogenizing.
    // The iterate stored internally also knows how to evaluate the gradient
    // efficiently, so our GradientEvaluator below just accesses it.
    struct DLibObjectiveEvaluator {
        DLibObjectiveEvaluator(_Inflator &inflator,
                WCStressOptimization::Objective<_N> &fullObjective)
            : m_iterate(NULL), m_inflator(inflator), m_fullObjective(fullObjective) {
            nParams = m_inflator.numParameters();
        };

        // Objective evaluation
        double operator()(const dlib_vector &x) const {
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);

            m_iterate = PatternOptimization::getIterate(m_iterate, m_inflator, nParams, &x_vec[0],
                                                        m_fullObjective);
            return m_iterate->evaluateJFull();
        }

        const _Iterate<_Sim> &currentIterate() const {
            assert(m_iterate); return *m_iterate;
        };

        size_t nParams;
    private:
        // Iterate is mutable so that operator() can be const as dlib requires
        mutable std::shared_ptr<_Iterate<_Sim>> m_iterate;
        _Inflator &m_inflator;
        WCStressOptimization::Objective<_N> &m_fullObjective;
    };
    
    // Extracts gradient from the iterate constructed by DLibObjectiveEvaluator.
    struct DLibGradientEvaluator {
        DLibGradientEvaluator(const DLibObjectiveEvaluator &obj) : m_obj(obj) { }

        // Gradient evaluation
        dlib_vector operator()(const dlib_vector &x) const {
            size_t nParams = m_obj.nParams;
            vector<Real> x_vec(nParams);
            for (size_t p = 0; p < nParams; ++p)
                x_vec[p] = x(p);
            if (m_obj.currentIterate().paramsDiffer(nParams, &x_vec[0]))
                throw std::runtime_error("Objective must be evaluated first");

            SField gradp = m_obj.currentIterate().gradp_JFull();

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
    void optimize_bfgs(SField &params, WCStressOptimization::Objective<_N> &fullObjective,
                       size_t niters, const string &outPath, size_t max_size) {
        DLibObjectiveEvaluator obj(m_inflator, fullObjective);
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
    _Inflator &m_inflator;
    std::vector<Real> m_patternParams, m_radiusBounds, m_transBounds, m_blendingBounds;
    std::map<size_t, Real> m_varLowerBounds, m_varUpperBounds;
};

}

#endif /* end of include guard: WCSOPTIMIZATION_HH */
