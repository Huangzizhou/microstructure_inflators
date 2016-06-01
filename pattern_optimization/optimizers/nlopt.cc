#include "nlopt.hh"

#include <nlopt.hpp>
#include <limits>

#include <vector>
#include <cassert>

using namespace PatternOptimization;

struct NLOptState {
    NLOptState(IterateManagerBase &im) : im(im), best(std::numeric_limits<Real>::max()) { }

    // TODO: how to choose when to report when there are constraints? We might
    // have to hack the NLOpt source.
    bool isImprovement(Real val) {
        if (val < best) {
            best = val;
            ++niters;
            return true;
        }
        return false;
    }

    size_t niters = 0;
    IterateManagerBase &im;
    Real best;
    std::string outPath;
};


double costFunc(const std::vector<double> &x, std::vector<double> &grad, void *optStateVoid)
{
    NLOptState *optState = reinterpret_cast<NLOptState *>(optStateVoid);
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
    if (optState->isImprovement(val)) {
        it.writeDescription(std::cout);
        std::cout << std::endl;
        if (optState->outPath != "")
            it.writeMeshAndFields(optState->outPath + "_" + std::to_string(optState->niters));
    }


    return it.evaluate();
}

void optimize_nlopt_sqslp(ScalarField<Real> &params,
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

    opt.set_maxeval(oconfig.niters);
    opt.set_xtol_rel(1e-16);

    // nlopt's C++ interface works on std::vectors
    std::vector<double> x(nParams);
    for (size_t p = 0; p < nParams; ++p)
        x[p] = params[p];

    double minf;
    /* nlopt::result result = */ opt.optimize(x, minf);

    // Convert the solution back.
    for (size_t p = 0; p < nParams; ++p)
        params[p] = x[p];
}

