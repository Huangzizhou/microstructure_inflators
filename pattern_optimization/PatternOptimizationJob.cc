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
//          "initial_params": [ ... ]
//      }
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  10/04/2014 17:15:02
////////////////////////////////////////////////////////////////////////////////
#include "PatternOptimizationJob.hh"
#include <stdexcept>
#include <fstream>

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

    if (dim == 2) {
        auto job2D = new Job<2>();
        job2D->targetMaterial.setFromPTree(materialSpec);
        parseVector(paramVals, job2D->initialParams);
        return job2D;
    }
    else if (dim == 3) {
        auto job3D = new Job<3>();
        job3D->targetMaterial.setFromPTree(materialSpec);
        parseVector(paramVals, job3D->initialParams);
        return job3D;
    }
    else throw runtime_error("Invalid dimension.");
}

}
