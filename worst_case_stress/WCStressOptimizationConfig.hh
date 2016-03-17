////////////////////////////////////////////////////////////////////////////////
// WCStressOptimizationConfig.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Simple configuration system hack (I don't want to pass these options
//      through all levels of the heirarchy).
//      This is a simple singleton pattern implementation that is NOT
//      threadsafe.
/      
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  12/16/2015 23:56:18
////////////////////////////////////////////////////////////////////////////////
#ifndef WCSTRESSOPTIMIZATIONCONFIG_HH
#define WCSTRESSOPTIMIZATIONCONFIG_HH

namespace WCStressOptimization {

class Config {
private:
    Config() { }
public:
    // By default, an "Lp norm" objective is really the p^th power of the Lp norm.
    // To use the true "Lp norm", globalObjectiveRoot must be set to
    // 2.0 * globalObjectivePNorm (since pointwise WCS is already squared (e.g. Frobenius) norm)
    Real globalObjectivePNorm = 1.0;
    Real globalObjectiveRoot  = 1.0;
    bool useVtxNormalPerturbationGradientVersion = false;
    static Config &get() {
        static Config configSingleton;
        return configSingleton;
    }

    // Whether to project out the normal stress component in the fluctuation
    // displacement shape derivative cell problem load.
    bool projectOutNormalStress = false;
};

}

#endif /* end of include guard: WCSTRESSOPTIMIZATIONCONFIG_HH */
