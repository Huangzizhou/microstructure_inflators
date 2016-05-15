#ifndef DLIB_HH
#define DLIB_HH

#include <Fields.hh>
#include <string>

#include "../IterateManagerBase.hh"
#include "../BoundConstraints.hh"
#include "../OptimizerConfig.hh"

// bfgs when oconfig.lbfgs_memory == 0, lbfgs otherwise.
void optimize_dlib_bfgs(ScalarField<Real> &params,
        const PatternOptimization::BoundConstraints &bds,
        PatternOptimization::IterateManagerBase &im,
        const PatternOptimization::OptimizerConfig &oconfig,
        const std::string &outPath);

#endif /* end of include guard: DLIB_HH */
