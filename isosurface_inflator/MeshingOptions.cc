#include "MeshingOptions.hh"

#include <iostream>
#include <set>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace std;
using boost::property_tree::ptree;

void MeshingOptions::load(const std::string &jsonPath) {
    ptree pt;
    read_json(jsonPath, pt);
    set<string> expectedKeys = {
        "domainErrorBound", "facetAngle", "facetSize", "facetDistance",
        "cellSize", "edgeSize", "cellRadiusEdgeRatio",
        "marchingSquaresGridSize", "marchingCubesGridSize",
        "maxArea", "featureAngleThreshold", "forceMSGridSize",
        "marchingSquaresCoarsening"
    };
    // Validate keys
    for (const auto &v : pt) {
        if (expectedKeys.count(v.first) == 0)
            cerr << "WARNING: ignoring unexpected meshing option " << v.first << endl;
    }

    // CGAL volume mesher options
    domainErrorBound        = pt.get<double>("domainErrorBound");
    facetAngle              = pt.get<double>("facetAngle");
    facetSize               = pt.get<double>("facetSize");
    facetDistance           = pt.get<double>("facetDistance");
    cellSize                = pt.get<double>("cellSize");
    edgeSize                = pt.get<double>("edgeSize");
    cellRadiusEdgeRatio     = pt.get<double>("cellRadiusEdgeRatio");
    marchingSquaresGridSize = pt.get<size_t>("marchingSquaresGridSize");
    marchingCubesGridSize   = pt.get<size_t>("marchingCubesGridSize");

    // 2D Mesher Options
    maxArea               = pt.get<double>("maxArea");
    featureAngleThreshold = pt.get<double>("featureAngleThreshold");

    // Optional flags
    if (pt.count("forceMSGridSize"))
        forceMSGridSize = pt.get<bool>("forceMSGridSize");
    if (pt.count("marchingSquaresCoarsening"))
        marchingSquaresCoarsening = pt.get<size_t>("marchingSquaresCoarsening");
}
