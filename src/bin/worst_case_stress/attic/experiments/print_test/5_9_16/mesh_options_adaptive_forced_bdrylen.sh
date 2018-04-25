maxArea=$1
msGridSize=$2
maxBdryEdgeLen=$3
cat <<HERE
{
    "domainErrorBound"            : 1e-5,
    "facetAngle"                  : 30.0,
    "facetSize"                   : 0.025,
    "facetDistance"               : 2e-3,
    "cellSize"                    : 0.15,
    "edgeSize"                    : 0.025,
    "cellRadiusEdgeRatio"         : 2.0,
    "marchingSquaresGridSize"     : $msGridSize,
    "marchingCubesGridSize"       : 128,

    "maxArea"                     : $maxArea,
    "featureAngleThreshold"       : 0.7853981633974483,
    "forceMSGridSize"             : true,
    "marchingSquaresCoarsening"   : 3,
    "curvatureAdaptive"           : true,
    "forceMaxBdryEdgeLen"         : $maxBdryEdgeLen
}
HERE