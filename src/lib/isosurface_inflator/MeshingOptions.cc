#include "MeshingOptions.hh"
#include <MeshFEM/StringUtils.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <set>

using json = nlohmann::json;

void MeshingOptions::load(const std::string &jsonPath) {
    std::ifstream is(jsonPath);
    if (!is.is_open()) {
        throw std::runtime_error("Cannot read json file: " + jsonPath);
    }
    json config;
    is >> config;
    load(config);

#if 0
    ptree pt;
    read_json(jsonPath, pt);
    std::set<std::string> expectedKeys = {
        "domainErrorBound", "facetAngle", "facetSize", "facetDistance",
        "cellSize", "edgeSize", "cellRadiusEdgeRatio",
        "marchingSquaresGridSize", "marchingCubesGridSize",
        "maxArea", "featureAngleThreshold", "forceMSGridSize",
        "marchingSquaresCoarsening", "curvatureAdaptive",
        "curveSimplifier",
        "forceMaxBdryEdgeLen",
        "jointBlendingMode",
        "jointBlendingFunction",
        "forceConsistentInterfaceMesh",
        "jacobian"
    };
    // Validate keys
    for (const auto &v : pt) {
        if (expectedKeys.count(v.first) == 0) {
            std::cerr << "WARNING: ignoring unexpected meshing option " << v.first << std::endl;
        }
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
        else if (boost::iequals(simp, "none"))
                curveSimplifier = NONE;
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

    if (pt.count("jointBlendingFunction")) {
        const std::string functionString = pt.get<std::string>("jointBlendingFunction");

        if (boost::iequals(functionString, "EXPONENTIAL")) {
            jointBlendingFunction = JointBlendFunction::EXPONENTIAL;
        }
        else if (boost::iequals(functionString, "POLY_SYMMETRIC")) {
            jointBlendingFunction = JointBlendFunction::POLY_SYMMETRIC;
        }
        else if (boost::iequals(functionString, "POLY_NONCONVEX")) {
            jointBlendingFunction = JointBlendFunction::POLY_NONCONVEX;
        }
        else if (boost::iequals(functionString, "POLY_PIECEWISE")) {
            jointBlendingFunction = JointBlendFunction::POLY_PIECEWISE;
        }
        else { throw std::runtime_error("Unrecognized blending function: " + functionString); }
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
#endif
}

void MeshingOptions::load(const nlohmann::json &config) {
    std::set<std::string> expectedKeys = {
        "domainErrorBound", "facetAngle", "facetSize", "facetDistance",
        "cellSize", "edgeSize", "cellRadiusEdgeRatio",
        "marchingSquaresGridSize", "marchingCubesGridSize",
        "maxArea", "featureAngleThreshold", "forceMSGridSize",
        "marchingSquaresCoarsening", "curvatureAdaptive",
        "curveSimplifier",
        "forceMaxBdryEdgeLen",
        "jointBlendingMode",
        "jointBlendingFunction",
        "forceConsistentInterfaceMesh",
        "jacobian"
    };

    // Validate keys
    for (auto it : config.items()) {
        if (expectedKeys.count(it.key()) == 0) {
            std::cerr << "WARNING: ignoring unexpected meshing option " << it.key() << std::endl;
        }
    }

    // CGAL volume mesher options
    domainErrorBound        = config["domainErrorBound"];
    facetAngle              = config["facetAngle"];
    facetSize               = config["facetSize"];
    facetDistance           = config["facetDistance"];
    cellSize                = config["cellSize"];
    edgeSize                = config["edgeSize"];
    cellRadiusEdgeRatio     = config["cellRadiusEdgeRatio"];
    marchingSquaresGridSize = config["marchingSquaresGridSize"];
    marchingCubesGridSize   = config["marchingCubesGridSize"];

    // 2D Mesher Options
    maxArea               = config["maxArea"];
    featureAngleThreshold = config["featureAngleThreshold"];

    // Optional flags
    forceMSGridSize = config.value("forceMSGridSize", forceMSGridSize);
    marchingSquaresCoarsening = config.value("marchingSquaresCoarsening", marchingSquaresCoarsening);
    curvatureAdaptive = config.value("curvatureAdaptive", curvatureAdaptive);
    if (config.count("forceMaxBdryEdgeLen")) {
        m_forceMaxBdryEdgeLen = true;
        m_forcedMaxBdryEdgeLen = config["forceMaxBdryEdgeLen"];
    }
    if (config.count("curveSimplifier")) {
        const std::string simp = MeshFEM::lowercase(config["curveSimplifier"]);
        if (simp == "collapse") {
            curveSimplifier = COLLAPSE;
        } else if (simp == "resample") {
            curveSimplifier = RESAMPLE;
        } else if (simp == "none") {
                curveSimplifier = NONE;
        } else {
            throw std::runtime_error("Unknown curve simplifier '" + simp + "'; expected 'collapse' or 'resample'");
        }
    }
    forceConsistentInterfaceMesh = config.value("forceConsistentInterfaceMesh", forceConsistentInterfaceMesh);

    if (config.count("jointBlendingMode")) {
        const std::string modeString = MeshFEM::lowercase(config["jointBlendingMode"]);
        if (modeString == "hull") {
            jointBlendingMode = JointBlendMode::HULL;
        } else if (modeString == "full") {
            jointBlendingMode = JointBlendMode::FULL;
        } else {
            throw std::runtime_error("Unrecognized blending mode: " + modeString);
        }
    }

    if (config.count("jointBlendingFunction")) {
        const std::string functionString = MeshFEM::lowercase(config["jointBlendingFunction"]);
        if (functionString == "exponential") {
            jointBlendingFunction = JointBlendFunction::EXPONENTIAL;
        }
        else if (functionString == "poly_symmetric") {
            jointBlendingFunction = JointBlendFunction::POLY_SYMMETRIC;
        }
        else if (functionString == "poly_nonconvex") {
            jointBlendingFunction = JointBlendFunction::POLY_NONCONVEX;
        }
        else if (functionString =="poly_piecewise") {
            jointBlendingFunction = JointBlendFunction::POLY_PIECEWISE;
        }
        else {
            throw std::runtime_error("Unrecognized blending function: " + functionString);
        }
    }

    if (config.count("jacobian")) {
        jacobian = read_jacobian(config["jacobian"]);
    } else {
        jacobian.setIdentity();
    }
}

Eigen::Matrix3d MeshingOptions::read_jacobian(const nlohmann::json &entry) {
    Eigen::Matrix3d jacobian = Eigen::Matrix3d::Identity();
    std::vector<double> x = entry;
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
    return jacobian;
}

