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
    Config() { }
public:
    bool ignoreShear = false;
    // For isosurface inflator: whether to use the normal velocity scalar field
    // directly (true) or to multiply by (vertex normal . face normal).
    // Note: the direct version empirically seems to be working better...
    // It should be used by default.
    bool useSDNormalShapeVelocityDirectly = true;
    ////////////////////////////////////////////////////////////////////////////
    // Remeshing configuration (for remeshing inflators)
    ////////////////////////////////////////////////////////////////////////////
    // Edge merging and splitting threshold (relative to median edge length).
    Real remeshMergeThreshold = 0.1;
    Real remeshSplitThreshold = 2.0;
    // 2D feature vertex detection threshold (on the angle difference from
    // M_PI, which corresponds to a perfectly straight boundary).
    Real remeshFeatureAngleThreshold = M_PI / 4;

    // Number of 1->4 triangle subdivision rounds before building
    // simulator. The subdivided geometry is used for the simulation mesh only;
    // everything else, like boundary parameter extraction, processes the
    // unsubdivided geometry.
    size_t fem2DSubdivRounds = 0;

    static Config &get() {
        static Config configSingleton;
        return configSingleton;
    }
};

}

#endif /* end of include guard: PATTERNOPTIMIZATIONCONFIG_HH */
