////////////////////////////////////////////////////////////////////////////////
// IterateFactory.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Create a PatternOptimization::Iterate (or sublass), and optionally
//  configure various objectives terms:
//      IFConfigTensorFit:               Fit to a specified compliance tensor.
//      IFConfigProximityRegularization: Regularize to stay close to init point.
//      IFConfigTargetVolume:            Target volume constraint.
//  (These particular configurations are provided by the ObjectiveTerm*.hh file)
//
//  This will require making IterateFactory and the individual term factories
//  templated by iterate type because templated virtual functions aren't
//  allowed.
//
//  Another option: make ObjectFactory<Iterate> inherit from a base class and
//  have getIterate() be a virtual function:
//      Problem: forwarding optional args to the iterate constructor is
//      impossible then. Maybe we don't need it if we have a single iterate
//      class, though.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  04/08/2016 00:15:28
////////////////////////////////////////////////////////////////////////////////
#ifndef ITERATEFACTORY_HH
#define ITERATEFACTORY_HH

#include "PatternOptimizationIterate.hh"
#include "ObjectiveTermNormalizations.hh"
#include "Inflator.hh"
#include <Future.hh>

namespace PatternOptimization {

// IterateFactory Configurations:
// These add objective terms to the iterate. The objectives can be configured
// by setting the corresponding public member variables inherited into the
// factory class.
// E.g.:
//  factory.IFConfigTargetVolume::enabled = true;
//  factory.IFConfigTargetVolume::targetVolume = 1.0;
//  factory.IFConfigTargetVolume::weight       = 1.0;
struct IFConfig {
    bool enabled = true;
};

////////////////////////////////////////////////////////////////////////////////
// IFApplyConfigs: Apply a list of IterateFactory configurations in sequence.
////////////////////////////////////////////////////////////////////////////////
template<class... IFConfigs>
struct IFApplyConfigs;

template<class _CF1, class... IFConfigs>
struct IFApplyConfigs<_CF1, IFConfigs...> {
    template<class _Factory, class _Iterate>
    static void apply(_Factory *factory, const std::unique_ptr<_Iterate> &it,
                      ObjectiveTermNormalizations &normalizations) {
        if (factory->_CF1::enabled) {
            _CF1 *config_ptr = factory;
            config_ptr->configIterate(it, normalizations);
        }
        IFApplyConfigs<IFConfigs...>::apply(factory, it, normalizations);
    }
};

template<>
struct IFApplyConfigs<> {
    template<class _Factory, class _Iterate>
    static void apply(_Factory *, const std::unique_ptr<_Iterate> &, ObjectiveTermNormalizations &) { }
};

template<class _Iterate, class... IFConfigs>
struct IterateFactory : public IFConfigs... {
    static constexpr size_t N = _Iterate::_N;
    using Iterate = _Iterate;
    using _Inflator = Inflator<N>;

    IterateFactory(_Inflator &inflator)
        : m_inflator(inflator) { }

    // Use previous iterate if evaluating the same point. Otherwise, attempt to
    // inflate the new parameters. Try three times to inflate, and if
    // unsuccessful, estimate the point with linear extrapolation.
    //
    // Because this is a variadic template, it can be used with an iterate
    // class whose constructor has arbitrary trailing arguments.
    std::unique_ptr<Iterate>
    getIterate(std::unique_ptr<Iterate> oldIterate, size_t nParams, const double *params) {
        if (oldIterate && !oldIterate->paramsDiffer(nParams, params)) {
            // Evaluating at same parameters as old iterate--old iterate is
            // exact and should already be configured/evaluated appropriately.
            oldIterate->disableEstimation();
            return oldIterate;
        }

        std::cerr << "INFLATING PARAMS:";
        for (size_t p = 0; p < nParams; ++p) { std::cerr << "\t" << params[p]; }
        std::cerr << std::endl;

        // Free up memory by releasing the old iterate if we aren't going to
        // use it for estimation
        if (!m_allowEstimation) { oldIterate.reset(); }

        std::unique_ptr<_Iterate> newIterate;
        bool success;
        for (size_t i = 0; i < m_numInflationAttempts; ++i) {
            success = true;
            try {
                newIterate = Future::make_unique<_Iterate>(m_inflator, nParams, params);
            }
            catch (std::exception &e) {
                std::cerr << "INFLATOR FAILED: " << e.what() << std::endl;
                success = false;
            }
            if (success) break;
        }
        if (!success) {
            std::cerr << "3 INFLATION ATTEMPTS FAILED." << std::endl;
            if (!m_allowEstimation) throw std::runtime_error("Inflation failure with estimation disabled");
            if (!oldIterate) throw std::runtime_error("Inflation failure on first iterate");
            // Extrapolate the old iterate to the new evaluation point.
            oldIterate->estimatePoint(nParams, params);
            return oldIterate;
        }

        // We actually created a new iterate--configure it accordingly.
        BENCHMARK_START_TIMER_SECTION("Adding objective terms");
        IFApplyConfigs<IFConfigs...>::apply(this, newIterate, m_normalizations);
        BENCHMARK_STOP_TIMER_SECTION("Adding objective terms");

        BENCHMARK_START_TIMER_SECTION("Caching objective term gradients");
        newIterate->evaluateObjectiveTerms(m_inflator);
        BENCHMARK_STOP_TIMER_SECTION("Caching objective term gradients");

        // // Write debug steepest descent field
        // {
        //     MSHFieldWriter debug("descent_debug.msh", newIterate->mesh());
        //     VectorField<Real, N> bdescent = m_inflator.boundaryVFieldFromParams(newIterate->steepestDescent());
        //     VectorField<Real, N> vdescent(newIterate->mesh().numVertices());
        //     vdescent.clear();
        //     for (auto v : newIterate->mesh().vertices()) {
        //         auto bv = v.boundaryVertex();
        //         if (!bv) continue;
        //         vdescent(v.index()) = bdescent(bv.index());
        //     }
        //     debug.addField("JFull Steepest Descent", vdescent);
        // }

        return newIterate;
    }

    size_t numParameters() const { return m_inflator.numParameters(); }
    bool    isParametric() const { return m_inflator.isParametric(); }

private:
    ObjectiveTermNormalizations m_normalizations;
    _Inflator &m_inflator;
    size_t m_numInflationAttempts = 1;
    // TODO: it should be possible to clear all factored matrices and still
    // support estimation...
    bool m_allowEstimation        = false;
};

// Inflator template parameter is last so that it can be inferred...
template<class _Iterate, class... IFConfigs>
std::unique_ptr<IterateFactory<_Iterate, IFConfigs...>> make_iterate_factory(Inflator<_Iterate::_N> &inflator) {
    return Future::make_unique<IterateFactory<_Iterate, IFConfigs...>>(inflator);
}

}

#endif /* end of include guard: ITERATEFACTORY_HH */
