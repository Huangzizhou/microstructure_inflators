#include "MeshingOptions.hh"

#include <iostream>
#include <set>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string.hpp>

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
        "marchingSquaresCoarsening", "curvatureAdaptive",
        "curveSimplifier",
        "forceMaxBdryEdgeLen",
        "jointBlendingMode",
        "forceConsistentInterfaceMesh",
        "jacobian"
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
    if (pt.count("curvatureAdaptive"))
        curvatureAdaptive = pt.get<bool>("curvatureAdaptive");
    if (pt.count("forceMaxBdryEdgeLen")) {
        m_forceMaxBdryEdgeLen = true;
        m_forcedMaxBdryEdgeLen = pt.get<double>("forceMaxBdryEdgeLen");
    }
    if (pt.count("curveSimplifier")) {
        auto simp = pt.get<std::string>("curveSimplifier");
        if (boost::iequals(simp, "collapse"))
            curveSimplifier = COLLAPSE;
        else if (boost::iequals(simp, "resample"))
            curveSimplifier = RESAMPLE;
        else throw std::runtime_error("Unknown curve simplifier '" + simp +
                                      "'; expected 'collapse' or 'resample'");
    }
    if (pt.count("forceConsistentInterfaceMesh")) {
        forceConsistentInterfaceMesh = pt.get<bool>("forceConsistentInterfaceMesh");
    }

    if (pt.count("jointBlendingMode")) {
        const std::string modeString = pt.get<std::string>("jointBlendingMode");

        if (boost::iequals(modeString, "HULL")) {
            jointBlendingMode = JointBlendMode::HULL;
        }
        else if (boost::iequals(modeString, "FULL")) {
            jointBlendingMode = JointBlendMode::FULL;
        }
        else { throw std::runtime_error("Unrecognized blending mode: " + modeString); }
    }
    if (pt.count("jacobian")) {
        std::vector<double> x;
        for (auto& item : pt.get_child("jacobian")) {
            x.push_back(item.second.get_value<double>());
        }
        if (x.size() == 4) {
            jacobian <<
                x[0], x[1], 0,
                x[2], x[3], 0,
                0   , 0   , 1;
        } else if (x.size() == 9) {
            // We read data as row-major matrix, but Eigen matrices default
            // to column-major storage, hence the transpose
            std::copy_n(x.data(), 9, jacobian.data());
            jacobian.transposeInPlace();
        } else {
            throw std::runtime_error("Invalid Jacobian matrix size: " + std::to_string(x.size()));
        }
    } else {
        jacobian.setIdentity();
    }
}
