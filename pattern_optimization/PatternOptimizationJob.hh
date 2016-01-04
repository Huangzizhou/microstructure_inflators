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
//          "initial_params": [ ... ],
//          "radiusBounds": [ low, high ],
//          "translationBounds": [ low, high ],
//          "paramConstraints": [ "p1 = p2 + 5", ... ],
//          "bounds": [{ "var": 0, "lower": 1, "upper": 2 }]
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
#include <map>
#include <string>
#include <stdexcept>
#include <fstream>
#include <boost/optional.hpp>

namespace PatternOptimization  {
class JobBase {
public:
    virtual ~JobBase() { }

    std::vector<Real> initialParams, radiusBounds, translationBounds;
    // The ground-truth parameters can be stored here--they are written to the
    // job file for reference.
    std::vector<Real> trueParams;
    std::vector<std::string> parameterConstraints;
    std::map<size_t, Real> varLowerBounds, varUpperBounds;
    boost::optional<Real> targetVolume;
};

template<size_t _N>
class Job : public JobBase {
public:
    size_t numParams() const { return initialParams.size(); }
    void writeJobFile(const std::string &jobFile) const {
        std::ofstream os(jobFile);
        if (!os.is_open()) 
            throw std::runtime_error("Couldn't open output job file " + jobFile);
        os << "{" << std::endl
           << "\t\"dim\": " << _N << "," << std::endl
           << "\t\"target\": " << targetMaterial << "," << std::endl;
        if (targetVolume)
            os << "\t\"target volume\": " << targetVolume << "," << std::endl;
        os << "\t\"initial_params\": [";
        for (size_t i = 0; i < initialParams.size(); ++i)
            os << (i ? ", " : "") << initialParams[i];
        os << "]," << std::endl;

        if (trueParams.size() == initialParams.size()) {
            os << "\t\"# true_params\": [";
            for (size_t i = 0; i < trueParams.size(); ++i)
                os << (i ? ", " : "") << trueParams[i];
            os << "]," << std::endl;
        }

        os << "\t\"radiusBounds\": [" << radiusBounds[0] << ", " << radiusBounds[1] << "]," << std::endl
           << "\t\"translationBounds\": [" << translationBounds[0] << ", " << translationBounds[1] << "]" << std::endl
           << "}" << std::endl;
    }

    virtual ~Job() { }

    Materials::Constant<_N> targetMaterial;
};

JobBase *parseJobFile(const std::string &jobFile);

}

#endif /* end of include guard: PATTERNOPTIMIZATIONJOB_HH */
