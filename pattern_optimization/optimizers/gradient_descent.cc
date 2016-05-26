#include "gradient_descent.hh"

void optimize_gd(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath)
{
    for (size_t i = 0; i < oconfig.niters; ++i) {
        auto &it = im.get(params.size(), params.data());
        it.writeDescription(std::cout);
        std::cout << std::endl;

        params += it.steepestDescent() * oconfig.gd_step;

        // Apply bound constraints
        for (size_t p = 0; p < im.numParameters(); ++p) {
            if (bds.hasUpperBound.at(p)) params[p] = std::min(params[p], bds.upperBound.at(p));
            if (bds.hasLowerBound.at(p)) params[p] = std::max(params[p], bds.lowerBound.at(p));
        }

        if (outPath != "") it.writeMeshAndFields(outPath + "_" + std::to_string(i));
    }
}
