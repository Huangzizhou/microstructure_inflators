#!/usr/bin/env python
import sys
num, = sys.argv[1:]

facetDistance = [1.6e-3, 1.1e-3, 8.0e-4, 5.6e-4, 4.0e-4, 2.8e-4, 2.0e-4,
                 1.6e-3, 1.1e-3, 8.0e-4, 5.6e-4, 4.0e-4, 2.8e-4, 2.0e-4]
edgeSize = [0.0320, 0.0226, 0.0160, 0.0113, 0.0080, 0.0056, 0.0040,
            0.0320, 0.0226, 0.0160, 0.0113, 0.0080, 0.0056, 0.0040]
cellSize = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
            0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]

i = int(num)

print('''
{
    "domainErrorBound"        : 1e-5,
    "facetAngle"              : 30.0,
    "facetSize"               : 0.02,
    "facetDistance"           : %f,
    "cellSize"                : %f,
    "edgeSize"                : %f,
    "cellRadiusEdgeRatio"     : 2.0,

    "marchingSquaresGridSize" : 2048,
    "marchingSquaresCoarsening": 3,
    "marchingCubesGridSize"   : 128,

    "maxArea"                 : 0.001,
    "featureAngleThreshold"   : 0.7853981633974483
}
''' % (facetDistance[i], cellSize[i], edgeSize[i]))
