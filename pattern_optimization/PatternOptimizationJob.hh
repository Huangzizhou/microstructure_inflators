////////////////////////////////////////////////////////////////////////////////
// PatternOptimizationJob.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Reads in JSON files specifying an optimization job. These job
//      specifications include a target material and the initial pattern
//      parameters from which to start optimization.
//      Format:
//      {
//          "dim": 2,
//          "target": {
//              <material_spec>
//          },
//          "initial_params": [ ... ]
//      }
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/04/2014 17:15:02
////////////////////////////////////////////////////////////////////////////////
#ifndef PATTERNOPTIMIZATIONJOB_HH
#define PATTERNOPTIMIZATIONJOB_HH
#include <Materials.hh>
#include <vector>
#include <string>
#include <stdexcept>

namespace PatternOptimization  {
class JobBase {
public:
    virtual ~JobBase() { }
};

template<size_t _N>
class Job : public JobBase {
public:
    size_t numParams() const { return initialParams.size(); }
    virtual ~Job() { }

    Materials::Constant<_N> targetMaterial;
    std::vector<Real> initialParams;
};

JobBase *parseJobFile(const std::string &jobFile);

}

#endif /* end of include guard: PATTERNOPTIMIZATIONJOB_HH */
