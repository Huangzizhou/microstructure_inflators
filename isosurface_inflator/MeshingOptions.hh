#ifndef MESHINGOPTIONS_HH
#define MESHINGOPTIONS_HH
#include <string>

struct MeshingOptions {
    MeshingOptions() { }
    MeshingOptions(const std::string &jsonPath) {
        load(jsonPath);
    }

    // CGAL volume mesher options
    double domainErrorBound          = 1e-5;
    double facetAngle                = 30.0;
    double facetSize                 = 0.025;
    double facetDistance             = 2e-3;
    double cellSize                  = 0.15;
    double edgeSize                  = 0.025;
    double cellRadiusEdgeRatio       = 2.0;

    // Cell face 2D mesher options
    size_t marchingSquaresGridSize = 256;

    // VCG Marching cubes options
    size_t marchingCubesGridSize = 128;

    void load(const std::string &jsonPath);

};

#endif /* end of include guard: MESHINGOPTIONS_HH */
