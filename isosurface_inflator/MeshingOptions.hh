#ifndef MESHINGOPTIONS_HH
#define MESHINGOPTIONS_HH
#include <string>
#include <cmath>

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

    // For MidplaneMesher, there is a single parameter: maxArea
    // The marchingSquaresGridSize and boundary max edge length are
    // automatically determined from this parameter:
    // The boundary max edge length is computed based on the edge length of an
    // equilateral triangle of area maxArea:
    //      maxLen = sqrt((4.0 / sqrt(3.0)) * maxArea)
    // Then, the marching squares grid size is chosen to ensure this maximal
    // edge length, even for diagonal edges. We aim for a square grid, so
    // we assume the diagonal is of length sqrt(2) * (side length) ==>
    //      sideLen = maxLen / sqrt(2)
    //      gridSize = ceil(domainLength / sideLen)
    double maxArea = 0.001;
    double featureAngleThreshold = M_PI / 4.0;
    // Derived parameters
    double maxEdgeLenFromMaxArea() const { return sqrt((4.0 / sqrt(3.0)) * maxArea); }
    size_t msGridSizeFromMaxArea(double domainLength) const {
        return ceil((domainLength * sqrt(2.0)) / maxEdgeLenFromMaxArea());
    }

    void load(const std::string &jsonPath);
};

#endif /* end of include guard: MESHINGOPTIONS_HH */
