#ifndef BOUNDCONSTRAINTS_HH
#define BOUNDCONSTRAINTS_HH

#include <vector>
#include <map>
#include "Inflator.hh"

namespace PatternOptimization {

struct BoundConstraints {
    // lb, ub: manually configured lower and upper bounds. All other bounds are
    // set using the defaults in radiusBounds[01] and translationBounds[01],
    // selected based on parameter type as determined by inflator.
    BoundConstraints(const InflatorBase &inflator,
                     const std::vector<Real> &radiusBounds,
                     const std::vector<Real> &translationBounds,
                     const std::vector<Real> &blendingBounds,
                     const std::map<size_t, Real> &lb,
                     const std::map<size_t, Real> &ub) {
        size_t nparams = inflator.numParameters();

        // All params get bounds by default.
        hasLowerBound.assign(nparams, true);
        hasUpperBound.assign(nparams, true);
        lowerBound.assign(nparams, 0.0);
        upperBound.assign(nparams, 0.0);

        // Explicitly specified bounds override default type-based bounds.
        // Set lower bounds
        for (size_t p = 0; p < nparams; ++p) {
            if (lb.count(p)) { lowerBound.at(p) = lb.at(p); continue; }
            switch (inflator.parameterType(p)) {
                case ParameterType::Thickness: lowerBound.at(p) =      radiusBounds.at(0); break;
                case ParameterType::Offset:    lowerBound.at(p) = translationBounds.at(0); break;
                case ParameterType::Blending:  lowerBound.at(p) =    blendingBounds.at(0); break;
                default: assert(false);
            }
        }

        for (size_t p = 0; p < nparams; ++p) {
            if (ub.count(p)) { upperBound.at(p) = ub.at(p); continue; }
            switch (inflator.parameterType(p)) {
                case ParameterType::Thickness: upperBound.at(p) =      radiusBounds.at(1); break;
                case ParameterType::Offset:    upperBound.at(p) = translationBounds.at(1); break;
                case ParameterType::Blending:  upperBound.at(p) =    blendingBounds.at(1); break;
                default: assert(false);
            }
        }
    }

    std::vector<Real> lowerBound, upperBound;
    std::vector<bool> hasLowerBound, hasUpperBound;
};

}

#endif /* end of include guard: BOUNDCONSTRAINTS_HH */
