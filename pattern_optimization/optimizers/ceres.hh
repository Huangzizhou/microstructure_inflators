#ifndef CERES_HH
#define CERES_HH

#include <Fields.hh>
#include <string>

#include "../IterateManagerBase.hh"
#include "../BoundConstraints.hh"
#include "../OptimizerConfig.hh"

void optimize_ceres_lm(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath);

// MHS on AUG 25, 2015:
// DOGLEG gets similar patterns (excluding rotations) for deformed cells with
// Jacobian=[a 0; 0 b] and Jacobian=[b 0; 0 a]
void optimize_ceres_dogleg(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath);

#endif /* end of include guard: CERES_HH */
