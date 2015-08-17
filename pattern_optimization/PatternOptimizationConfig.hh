////////////////////////////////////////////////////////////////////////////////
// PatternOptimizationConfig.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Simple configuration system hack (I don't want to pass these options
//      through all levels of the heirarchy).
//      This is a simple singleton pattern implementation that is NOT
//      threadsafe.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  08/11/2015 18:30:35
////////////////////////////////////////////////////////////////////////////////
#ifndef PATTERNOPTIMIZATIONCONFIG_HH
#define PATTERNOPTIMIZATIONCONFIG_HH

namespace PatternOptimization {

class Config {
private:
    Config() : ignoreShear(false) { }
public:
    bool ignoreShear;
    static Config &get() {
        static Config configSingleton;
        return configSingleton;
    }
};

}

#endif /* end of include guard: PATTERNOPTIMIZATIONCONFIG_HH */
