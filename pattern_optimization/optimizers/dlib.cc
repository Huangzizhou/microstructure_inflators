#include "dlib.hh"

#include <dlib/optimization.h>
#include <vector>

using dlib_vector = dlib::matrix<double, 0, 1>;
using namespace PatternOptimization;

// Evaluates the objective by inflating the wire mesh and homogenizing.
// The iterate stored internally also knows how to evaluate the gradient
// efficiently, so our GradientEvaluator below just accesses it.
struct DLibObjectiveEvaluator {
    DLibObjectiveEvaluator(IterateManagerBase &imanager) : m_imanager(imanager) { };

    double operator()(const dlib_vector &x) const {
        std::vector<Real> x_vec(m_imanager.numParameters());
        for (size_t p = 0; p < x_vec.size(); ++p)
            x_vec[p] = x(p);

        auto &it = m_imanager.get(x_vec.size(), &x_vec[0]);
        return it.evaluate();
    }
private:
    // Iterate manager is mutable so that operator() can be const as dlib requires
    IterateManagerBase &m_imanager;
};

// Extracts gradient from the iterate constructed by DLibObjectiveEvaluator.
struct DLibGradientEvaluator {
    DLibGradientEvaluator(const IterateManagerBase &imanager) : m_imanager(imanager) { }

    dlib_vector operator()(const dlib_vector &x) const {
        std::vector<Real> x_vec(m_imanager.numParameters());
        for (size_t p = 0; p < x_vec.size(); ++p)
            x_vec[p] = x(p);
        if (m_imanager.get().paramsDiffer(x_vec.size(), &x_vec[0]))
            throw std::runtime_error("Objective must be evaluated first");
        ScalarField<Real> gradp(m_imanager.get().gradp());

        dlib_vector result(x_vec.size());
        for (size_t p = 0; p < x_vec.size(); ++p)
            result(p) = gradp[p];
        return result;
    }
private:
    const IterateManagerBase &m_imanager;
};

// Hack to get notified at the end of each iteration: subclass the stop
// strategy.
class ReportingStopStrategy : public dlib::objective_delta_stop_strategy {
    typedef dlib::objective_delta_stop_strategy Base;
public:
    ReportingStopStrategy(double min_delta, unsigned long max_iter,
                          const IterateManagerBase &imanager, const std::string &outPath)
        : Base(min_delta, max_iter), m_imanager(imanager), m_outPath(outPath), m_iter(0) { }

    template <typename T>
    bool should_continue_search(const T& x, const double funct_value,
        const T& funct_derivative) {
        m_imanager.get().writeDescription(std::cout);
        std::cout << std::endl;
        if (m_outPath != "")
            m_imanager.get().writeMeshAndFields(m_outPath + "_" + std::to_string(m_iter));
        ++m_iter;
        return Base::should_continue_search(x, funct_value, funct_derivative);
    }

private:
    const IterateManagerBase &m_imanager;
    std::string m_outPath;
    size_t m_iter;
};

void optimize_dlib_bfgs(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath)
{
    if (!im.isParametric()) throw std::runtime_error("BFGS only works with parametric inflators");
    DLibObjectiveEvaluator obj(im);
    DLibGradientEvaluator grad(im);

    size_t nParams = im.numParameters();
    // convert initial parameter vector
    dlib_vector optParams(nParams);
    for (size_t p = 0; p < nParams; ++p)
        optParams(p) = params[p];

    dlib_vector lowerBounds(nParams), upperBounds(nParams);
    for (size_t p = 0; p < nParams; ++p) {
        if (!(bds.hasLowerBound.at(p) && bds.hasUpperBound.at(p)))
            throw std::runtime_error("Currently all parameters must be bounded for bfgs.");
        lowerBounds(p) = bds.lowerBound.at(p);
        upperBounds(p) = bds.upperBound.at(p);
    }

    size_t max_size = oconfig.lbfgs_memory;
    if (max_size == 0) {
        dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
                ReportingStopStrategy(1e-16, oconfig.niters, im, outPath),
                obj, grad, optParams, lowerBounds, upperBounds);
    }
    else {
        dlib::find_min_box_constrained(dlib::lbfgs_search_strategy(max_size),
                ReportingStopStrategy(1e-16, oconfig.niters, im, outPath),
                obj, grad, optParams, lowerBounds, upperBounds);
    }

    // convert solution
    for (size_t p = 0; p < nParams; ++p)
        params[p] = optParams(p);
}
