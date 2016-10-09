#include "nlopt.hh"

#include <nlopt.hpp>
#include <limits>

#include <vector>
#include <cassert>

#include "../pattern_optimization/Constraint.hh"

using namespace PatternOptimization;

struct NLOptState {
    NLOptState(IterateManagerBase &im) : im(im), best(std::numeric_limits<Real>::max()) { }

    // TODO: When constraints are involved, the iterate can be an improvement if
    // either infeasibility decreases or objective decreases.
    // (Currently just always reports)
    bool isImprovement(Real /* val */, const std::vector<Real> &params) {
        bool differ = true;
        if (prevParams.size() == params.size()) {
            differ = false;
            for (size_t i = 0; i < params.size(); ++i) {
                if (std::abs(prevParams[i] - params[i]) > 1e-9) {
                    differ = true;
                    break;
                }
            }
        }
        prevParams = params;
        if (!differ) return false;
#if 0
        if (val < best) {
            best = val;
            ++niters;
            return true;
        }
        return false;
#endif
		++niters;
        return true;
    }

    std::vector<Real> prevParams;

    size_t niters = 0;
    IterateManagerBase &im;
    Real best;
    std::string outPath;
};

struct NLOptConstraintEvaluator {
    NLOptConstraintEvaluator(NLOptState &s, size_t idx)
        : state(s), constraintIndex(idx) { }

    NLOptState &state;
    size_t constraintIndex;
};

double costFunc(const std::vector<double> &x, std::vector<double> &grad, void *optStateVoid) {
    auto optState = reinterpret_cast<NLOptState *>(optStateVoid);
    assert(optState);

    auto &it = optState->im.get(x.size(), &x[0]);

    if (!grad.empty()) {
        assert(grad.size() == x.size());
        ScalarField<Real> gp = it.gradp();
        assert(gp.domainSize() == grad.size());
        for (size_t p = 0; p < gp.domainSize(); ++p)
            grad[p] = gp[p];
    }

    Real val = it.evaluate();
    if (optState->isImprovement(val, x)) {
        it.writeDescription(std::cout);
        std::cout << std::endl;
        if (optState->outPath != "")
            it.writeMeshAndFields(optState->outPath + "_" + std::to_string(optState->niters));
    }


    return it.evaluate();
}

void constraintFunc(unsigned m, double *result, unsigned n, const double *x,
                    double *gradient, void *cevalVoid) {
    auto ceval = reinterpret_cast<NLOptConstraintEvaluator *>(cevalVoid);
    assert(ceval);

    auto &it = ceval->state.im.get(n, x);
    const EvaluatedConstraint &c = it.evaluatedConstraint(ceval->constraintIndex);
    assert(c.values.domainSize() == m);

    for (size_t i = 0; i < m; ++i)
        result[i] = c.values[i];

    // Unless the gradient is requested, we are done.
    if (gradient == NULL) return;

    for (size_t i = 0; i < m; ++i)
        for (size_t j = 0; j < n; ++j)
            gradient[i * n + j] = c.jacobian(i, j);
}

void optimize_nlopt_slsqp(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath)
{
    nlopt::opt opt(nlopt::LD_SLSQP, params.domainSize());

    const size_t nParams = params.domainSize();
    for (size_t p = 0; p < nParams; ++p) {
        if (!(bds.hasLowerBound.at(p) && bds.hasUpperBound.at(p)))
            throw std::runtime_error("Currently all parameters must be bounded.");
    }

    opt.set_lower_bounds(bds.lowerBound);
    opt.set_upper_bounds(bds.upperBound);

    NLOptState state(im);
    state.outPath = outPath;
    opt.set_min_objective(costFunc, (void *) &state);

    // Must create iterate to query the constraints.
    // This iterate for the initial parameters will be reused by the first
    // optimization iteration.
    std::vector<std::unique_ptr<NLOptConstraintEvaluator>> cevals;
    auto &it = im.get(nParams, &params[0]);
    for (size_t i = 0; i < it.numConstraints(); ++i) {
        const auto &c = it.evaluatedConstraint(i);
        const size_t m = c.dimension();
        // TODO: change to 1e-2
        std::vector<Real> tol(m, 1e-4);
        cevals.push_back(Future::make_unique<NLOptConstraintEvaluator>(state, i));
        auto ceval = cevals.back().get();
        switch (c.type) {
            case ConstraintType::  EQUALITY: opt.  add_equality_mconstraint(constraintFunc, (void *) ceval, tol); break;
            case ConstraintType::INEQUALITY: opt.add_inequality_mconstraint(constraintFunc, (void *) ceval, tol); break;
            default: throw std::runtime_error("Unexpected constraint type");
        }
    }

    opt.set_maxeval(oconfig.niters);
    opt.set_xtol_rel(1e-16);

    // nlopt's C++ interface works on std::vectors
    std::vector<double> x(nParams);
    for (size_t p = 0; p < nParams; ++p)
        x[p] = params[p];

    double minf = 0;
    /* nlopt::result result = */ opt.optimize(x, minf);

    // Convert the solution back.
    for (size_t p = 0; p < nParams; ++p)
        params[p] = x[p];
}

