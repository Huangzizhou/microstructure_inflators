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
    Real globalObjectivePNorm = 1.0;
    static Config &get() {
        static Config configSingleton;
        return configSingleton;
    }
};

}

#endif /* end of include guard: WCSTRESSOPTIMIZATIONCONFIG_HH */
