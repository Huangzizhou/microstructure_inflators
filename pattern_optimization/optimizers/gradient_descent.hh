#ifndef GRADIENT_DESCENT_HH
#define GRADIENT_DESCENT_HH

#include <Fields.hh>
#include <string>

#include "../IterateManagerBase.hh"
#include "../BoundConstraints.hh"
#include "../OptimizerConfig.hh"

void optimize_gd(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath);


#endif /* end of include guard: GRADIENT_DESCENT_HH */
