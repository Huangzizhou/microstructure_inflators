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

#include "Inflator.hh"

namespace PatternOptimization  {
class JobBase {
public:
    virtual ~JobBase() { }

    size_t numParams() const { return initialParams.size(); }

    // Verifies that the correct number of parameters were specified in the job
    // (must match inflator). For non-parametric inflators (like the
    // BoundaryPerturbationInflator) an incorrect number of parameters is
    // tolerated with a warning, in which case the initial parameters are taken
    // to be all zero.
    // Returns the (possibly modified) initial parameters
    std::vector<Real> validatedInitialParams(const InflatorBase &inflator) const {
        std::vector<Real> params(initialParams);
        // Allow non-parametric inflator to ignore initialParams (if size mismatch)
        if (!inflator.isParametric()) {
            if (numParams() != inflator.numParameters()) {
                if (numParams() > 0)
                    std::cerr << "WARNING: ignoring incorrectly-sized initial parameters for non-parametric inflator." << std::endl;
                params.assign(inflator.numParameters(), 0.0);
            }
        }
        if (numParams() != inflator.numParameters()) {
            for (size_t i = 0; i < inflator.numParameters(); ++i)
                std::cerr << "param " << i
                          << " role: " << parameterTypeString(inflator.parameterType(i))
                          << std::endl;
            throw std::runtime_error("Invalid number of parameters.");
        }

        for (const auto &boundEntry : varLowerBounds) {
            if (boundEntry.first > params.size())
                std::cerr << "WARNING: bound on nonexistent variable" << std::endl;
        }

        for (size_t p = 0; p < params.size(); ++p) {
            if (varLowerBounds.count(p)) {
                 if ((params[p] < varLowerBounds.at(p)) ||
                     (params[p] > varUpperBounds.at(p))) {
                    throw std::runtime_error("Initial point infeasible");
                 }
            }
        }
        return params;
    }


    std::vector<Real> initialParams, radiusBounds, translationBounds;
    std::vector<Real> blendingBounds = { 10.0, 100.0 };
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
    void writeJobFile(const std::string &jobFile) const {
        std::ofstream os(jobFile);
        if (!os.is_open()) 
            throw std::runtime_error("Couldn't open output job file " + jobFile);
        os << "{" << std::endl
           << "\t\"dim\": " << _N << "," << std::endl
           << "\t\"target\": " << targetMaterial << "," << std::endl;
        if (targetVolume)
            os << "\t\"target volume\": " << *targetVolume << "," << std::endl;
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
           << "\t\"blendingBounds\": [" << blendingBounds[0] << ", " << blendingBounds[1] << "]" << std::endl
           << "}" << std::endl;
    }

    virtual ~Job() { }

    Materials::Constant<_N> targetMaterial;
};

JobBase *parseJobFile(const std::string &jobFile);

}

#endif /* end of include guard: PATTERNOPTIMIZATIONJOB_HH */
