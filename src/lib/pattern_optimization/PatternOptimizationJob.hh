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

#include <MeshFEM/Materials.hh>
#include <inflators/Inflator.hh>
#include <boost/optional.hpp>
#include <vector>
#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include <fstream>
#include <iomanip>

namespace PatternOptimization  {
class JobBase {
public:
    virtual ~JobBase() { }

    size_t numParams() const {
        if (paramsMask.empty()) {
            return initialParams.size();
        }
        else {
            assert(paramsMask.size() == initialParams.size());
            return std::count(paramsMask.begin(), paramsMask.begin()+paramsMask.size(), false);
        }
    }

    // Verifies that the correct number of parameters were specified in the job
    // (must match inflator). For non-parametric inflators (like the
    // BoundaryPerturbationInflator) an incorrect number of parameters is
    // tolerated with a warning, in which case the initial parameters are taken
    // to be all zero.
    // Returns the (possibly modified) initial parameters
    std::vector<Real> validatedInitialParams(const InflatorBase &inflator) const {
        std::vector<Real> params(initialParams);
        // Set params to the default if they're omitted in the job file
        size_t np;
        if (numParams() == 0) {
            params = inflator.defaultParameters();
            np = params.size();
        }
        else {
            np = numParams();
        }
        //const size_t np = params.size();
        // Allow non-parametric inflator to ignore initialParams (if size mismatch)
        if (!inflator.isParametric()) {
            if (np != inflator.numParameters()) {
                std::cerr << "WARNING: ignoring incorrectly-sized initial parameters for non-parametric inflator." << std::endl;
                params.assign(inflator.numParameters(), 0.0);
            }
            return params;
        }
        if (np != inflator.numParameters()) {
            for (size_t i = 0; i < inflator.numParameters(); ++i)
                std::cerr << "param " << i
                          << " role: " << parameterTypeString(inflator.parameterType(i))
                          << std::endl;
            std::cerr <<  "Inflator was expecting " << inflator.numParameters() << " parameters, but input contained only "
                      << np << " values" << std::endl;
            throw std::runtime_error("Invalid number of parameters.");
        }

        for (const auto &boundEntry : varLowerBounds) {
            if (boundEntry.first > params.size())
                std::cerr << "WARNING: bound on nonexistent variable" << std::endl;
        }

        for (size_t p = 0; p < params.size(); ++p) {
            if (varLowerBounds.count(p)) {
                if (params[p] < varLowerBounds.at(p)) {
                    std::cerr << "WARNING: param " << p << " clamped to lower bound." << std::endl;
                    params[p] = varLowerBounds.at(p);
                }
                if (params[p] > varUpperBounds.at(p)) {
                    std::cerr << "WARNING: param " << p << " clamped to upper bound." << std::endl;
                    params[p] = varUpperBounds.at(p);
                }
            }
        }

        return params;
    }

    void writeJobFile(const std::string &jobFile) const {
        std::ofstream os(jobFile);
        if (!os.is_open())
            throw std::runtime_error("Couldn't open output job file " + jobFile);
        writeJobFile(os);
    }

    virtual void writeJobFile(std::ostream &os) const = 0;

    std::vector<Real> initialParams, radiusBounds, translationBounds;
    std::vector<bool> paramsMask;
    std::vector<std::string> metaParams;
    std::vector<Real> blendingBounds = { 10.0, 100.0 };
    std::vector<Real> metaBounds = { 0.01, 0.99 };
    std::vector<Real> custom1Bounds = { 0.01, 0.99 };
    std::vector<Real> custom2Bounds = { 0.01, 0.99 };
    std::vector<Real> custom3Bounds = { 0.01, 0.99 };
    std::vector<Real> custom4Bounds = { 0.01, 0.99 };
    std::vector<Real> custom5Bounds = { 0.01, 0.99 };
    std::vector<Real> custom6Bounds = { 0.01, 0.99 };
    std::vector<Real> custom7Bounds = { 0.01, 0.99 };
    std::vector<Real> custom8Bounds = { 0.01, 0.99 };
    // The ground-truth parameters can be stored here--they are written to the
    // job file for reference.
    std::vector<Real> trueParams;
    std::vector<std::string> parameterConstraints;
    std::map<size_t, Real> varLowerBounds, varUpperBounds;
    boost::optional<Real> targetVolume;
    size_t numberCustomTypes = 0;
};

template<size_t _N>
class Job : public JobBase {
public:
    virtual void writeJobFile(std::ostream &os) const {
        os << std::setprecision(19);
        os << "{" << std::endl
           << "\t\"dim\": " << _N << "," << std::endl
           << "\t\"target\": " << targetMaterial << "," << std::endl;
        if (targetVolume)
            os << "\t\"targetVolume\": " << *targetVolume << "," << std::endl;
        if (initialParams.size()) {
            os << "\t\"initial_params\": [";
            for (size_t i = 0; i < initialParams.size(); ++i)
                os << (i ? ", " : "") << std::setprecision(10) << initialParams[i];
            os << "]," << std::endl;
        }
        if (paramsMask.size()) {
            os << "\t\"paramsMask\": [";
            for (size_t i = 0; i < paramsMask.size(); ++i)
                os << (i ? ", " : "") << paramsMask[i];
            os << "]," << std::endl;
        }

        if (trueParams.size() == initialParams.size()) {
            os << "\t\"# true_params\": [";
            for (size_t i = 0; i < trueParams.size(); ++i)
                os << (i ? ", " : "") << trueParams[i];
            os << "]," << std::endl;
        }

        if (parameterConstraints.size()) {
            os << "\t\"paramConstraints\": [";
            bool first = true;
            for (const std::string &pc : parameterConstraints) {
                if (!first) os << ",";
                first = false;
                os << std::endl << "\t\t\"" << pc << "\"";
            }
            os << "\t]," << std::endl;
        }

        if (varLowerBounds.size() + varUpperBounds.size()) {
            os << "\t\"bounds\": [";
            bool first = true;
            for (size_t p = 0; p < initialParams.size(); ++p) {
                if (varLowerBounds.count(p) + varUpperBounds.count(p) == 0)
                    continue;
                if (!first) os << ",";
                first = false;
                os << std::endl << "\t\t{\"var\": " << p;
                if (varLowerBounds.count(p)) os << ", \"lower\": " << varLowerBounds.at(p);
                if (varUpperBounds.count(p)) os << ", \"upper\": " << varUpperBounds.at(p);
                os << "}";
            }
            os << std::endl << "\t]," << std::endl;
        }

        if (numberCustomTypes > 0)
            os << "\t\"custom1Bounds\": [" << custom1Bounds[0] << ", " << custom1Bounds[1] << "]," << std::endl;
        if (numberCustomTypes > 1)
            os << "\t\"custom2Bounds\": [" << custom2Bounds[0] << ", " << custom2Bounds[1] << "]," << std::endl;
        if (numberCustomTypes > 2)
            os << "\t\"custom3Bounds\": [" << custom3Bounds[0] << ", " << custom3Bounds[1] << "]," << std::endl;
        if (numberCustomTypes > 3)
            os << "\t\"custom4Bounds\": [" << custom4Bounds[0] << ", " << custom4Bounds[1] << "]," << std::endl;
        if (numberCustomTypes > 4)
            os << "\t\"custom5Bounds\": [" << custom5Bounds[0] << ", " << custom5Bounds[1] << "]," << std::endl;
        if (numberCustomTypes > 5)
            os << "\t\"custom6Bounds\": [" << custom6Bounds[0] << ", " << custom6Bounds[1] << "]," << std::endl;
        if (numberCustomTypes > 6)
            os << "\t\"custom7Bounds\": [" << custom7Bounds[0] << ", " << custom7Bounds[1] << "]," << std::endl;
        if (numberCustomTypes > 7)
            os << "\t\"custom8Bounds\": [" << custom8Bounds[0] << ", " << custom8Bounds[1] << "]," << std::endl;

        os << "\t\"radiusBounds\": [" << radiusBounds[0] << ", " << radiusBounds[1] << "]," << std::endl
           << "\t\"translationBounds\": [" << translationBounds[0] << ", " << translationBounds[1] << "]," << std::endl
           << "\t\"blendingBounds\": [" << blendingBounds[0] << ", " << blendingBounds[1] << "]" << std::endl;

        os << "}" << std::endl;
    }

    virtual ~Job() { }

    Materials::Constant<_N> targetMaterial;
};

std::unique_ptr<JobBase> parseJobFile(const std::string &jobFile);

}

#endif /* end of include guard: PATTERNOPTIMIZATIONJOB_HH */
