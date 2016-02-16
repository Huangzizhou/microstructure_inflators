////////////////////////////////////////////////////////////////////////////////
// PatternOptimizationJob.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Reads in JSON files specifying an optimization job. These job
//      specifications include a target material and the initial pattern
//      parameters from which to start optimization.
//      {
//          "dim": 2,
//          "target": {
//              <material_spec>
//          },
//          "initial_params": [ ... ],
//          "radiusBounds": [ low, high ],
//          "translationBounds": [ low, high ],
//          "paramConstraints": [ "p1 = p2 + 5", ... ],
//          "bounds": [{ "var": 0, "lower": 1, "upper": 2 }]
//      }
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/04/2014 17:15:02
////////////////////////////////////////////////////////////////////////////////
#include "PatternOptimizationJob.hh"
#include <stdexcept>
#include <fstream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace std;
using boost::property_tree::ptree;

namespace PatternOptimization {

void parseVector(const ptree &pt, vector<Real> &v) {
    v.clear();
    for (const auto &val : pt) {
        if (!val.first.empty()) throw runtime_error("Failed to parse vector");
        v.push_back(val.second.get_value<Real>());
    }
}

JobBase *parseJobFile(const string &jobFile) {
    ifstream is(jobFile);
    if (!is.is_open())
        throw runtime_error("Couldn't open job file " + jobFile);
    ptree pt;
    read_json(is, pt);

    size_t dim = pt.get<size_t>("dim");
    auto materialSpec = pt.get_child("target");
    auto paramVals = pt.get_child("initial_params");
    auto radiusBounds = pt.get_child("radiusBounds");
    auto translationBounds = pt.get_child("translationBounds");

    JobBase *job;
    if (dim == 2) {
        auto job2D = new Job<2>();
        job2D->targetMaterial.setFromPTree(materialSpec);
        job = job2D;
    }
    else if (dim == 3) {
        auto job3D = new Job<3>();
        job3D->targetMaterial.setFromPTree(materialSpec);
        job = job3D;
    }
    else throw runtime_error("Invalid dimension.");

    parseVector(paramVals, job->initialParams);
    parseVector(radiusBounds, job->radiusBounds);
    parseVector(translationBounds, job->translationBounds);

    if (pt.count("blendingBounds"))
        parseVector(pt.get_child("blendingBounds"), job->blendingBounds);

    if (pt.count("paramConstraints") != 0) {
        auto constraints = pt.get_child("paramConstraints");
        for (const auto &val : constraints) {
            if (!val.first.empty()) throw runtime_error("Failed to read constraints");
            job->parameterConstraints.emplace_back(val.second.get_value<string>());
        }
    }

    // individual variable overriding radius/translation bounds
    if (pt.count("bounds") !=  0) {
        try {
            auto bounds = pt.get_child("bounds");
            for (const auto &bound : bounds) {
                size_t var = bound.second.get<size_t>("var");
                job->varLowerBounds[var] = bound.second.get<Real>("lower");
                job->varUpperBounds[var] = bound.second.get<Real>("upper");
            }
        }
        catch (...) {
            throw std::runtime_error("Couldn't parse variable bounds.");
        }
    }

    job->targetVolume = pt.get_optional<double>("target volume");

    return job;
}

}
