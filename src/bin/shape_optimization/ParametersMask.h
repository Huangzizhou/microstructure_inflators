////////////////////////////////////////////////////////////////////////////////
// ParametersMask.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Given a pattern, a set of parameters and a set of boundary conditions,
//    decide which variables should not change during optimization.
//
*/
//  Author:  Davi Colli Tozoni (dctozoni) davi.tozoni@nyu.edu
//  Company:  New York University
//  Created:  2/13/18
////////////////////////////////////////////////////////////////////////////////

#ifndef MICROSTRUCTURES_PARAMETERSMASK_H
#define MICROSTRUCTURES_PARAMETERSMASK_H

#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <MeshFEM/Future.hh>
#include <isosurface_inflator/WireMesh.hh>
#include <inflators/wrappers/BoundaryPerturbationInflator.hh>
#include <isosurface_inflator/Symmetry.hh>

// Json
#include "json.hpp"
using json = nlohmann::json;

using namespace std;

namespace po = boost::program_options;

using Point = Point3<double>;

namespace ParametersMask {
     void jsonToRegions(std::string jsonPath, const BBox<Point> &bbox, vector<Region<Point> *> &inRegions, vector<Region<Point> *> &exceptRegions) {
        // reading JSON file
        ifstream input(jsonPath);
        json j;
        input >> j;

        json regions = j["regions"];
        // iterate the array
        for (json::iterator it = regions.begin(); it != regions.end(); ++it) {
            Point minCorner;
            Point maxCorner;

            json current_region = *it;

            json box = current_region["box"];
            json box_percentage = current_region["box%"];

            if (box.size() > 0) {
                vector<float> min_corner = box["minCorner"];
                vector<float> max_corner = box["maxCorner"];

                minCorner << min_corner[0], min_corner[1], min_corner[2];
                maxCorner << max_corner[0], max_corner[1], max_corner[2];
            }
            if (box_percentage.size() > 0) {
                // Find the bounding box
                Eigen::Vector3d m = bbox.minCorner;
                Eigen::Vector3d M = bbox.maxCorner;

                vector<float> min_corner = box_percentage["minCorner"];
                vector<float> max_corner = box_percentage["maxCorner"];

                Eigen::Array3d min_corner_array;
                Eigen::Array3d max_corner_array;
                min_corner_array << min_corner[0], min_corner[1], min_corner[2];
                max_corner_array << max_corner[0], max_corner[1], max_corner[2];

                // scale (w/2, h/2)
                double w = M(0) - m(0);
                double h = M(1) - m(1);
                double d = M(2) - m(2);
                min_corner_array(0) *= w;///2;
                min_corner_array(1) *= h;///2;
                min_corner_array(2) *= d;///2;
                max_corner_array(0) *= w;///2;
                max_corner_array(1) *= h;///2;
                max_corner_array(2) *= d;///2;

                // translate to min corner
                min_corner_array += m.array();
                max_corner_array += m.array();

                minCorner = min_corner_array.matrix();
                maxCorner = max_corner_array.matrix();
            }

            string typeString = current_region["type"];
            if (typeString.compare("dirichlet") == 0 || typeString.compare("force") == 0 || typeString.compare("zero") == 0) {
                Region<Point> * newRegion = new BBox<Point>(minCorner, maxCorner);
                exceptRegions.push_back(newRegion);
            }
            else if (typeString.compare("optimization") == 0) {
                Region<Point> * newRegion = new BBox<Point>(minCorner, maxCorner);
                inRegions.push_back(newRegion);
            }
            else {
                throw std::runtime_error("Region of type " + typeString + "not expected.");
            }
        }
    }

    vector<Point> excludedPoints(vector<Point> points, vector<Region<Point> *> &inRegions, vector<Region<Point> *> &exceptRegions) {
        set<int> in;

        // marking which points are filtered and which are not
        if (inRegions.size() > 0) {
            for (auto region : inRegions) {
                for (size_t i = 0; i < points.size(); i++) {
                    if (region->containsPoint(points[i])) {
                        in.insert(i);
                    }
                }
            }
        }
        else {
            for (size_t i = 0; i < points.size(); i++) {
                in.insert(i);
            }
        }

        for (auto region : exceptRegions) {
            for (size_t i = 0; i < points.size(); i++) {
                if (region->containsPoint(points[i])) {
                    in.erase(i);
                }
            }
        }

        vector<Point> result;
        for (size_t i = 0; i < points.size(); i++) {
            if (!in.count(i))
                result.push_back(points[i]);
        }

        return result;
    }

    void extractFilteringRegions(string bcondsPath, vector<Point> points, vector<Region<Point> *> &inRegions, vector<Region<Point> *> &exceptRegions) {
        if (!bcondsPath.empty()) {
            BBox<Point> bb(points);

            jsonToRegions(bcondsPath, bb, inRegions, exceptRegions);
        }
    }

    vector<bool> generateParametersMask(string patternPath, vector<double> params, string bcondsPath, size_t blendingPolySize = 0, string sym = "non_periodic") {
        WireMeshBase * wireMesh;

        vector<Point3<double>> points;
        vector<pair<size_t, size_t>> edges;
        vector<double> thicknesses;
        vector<double> blendingParams;
        vector<vector<double>> blendingPolyParams;

        if      (sym == "cubic"          ) {
            wireMesh = new WireMesh<Symmetry::Cubic<>>(patternPath, 0, blendingPolySize);
        }
        else if (sym == "orthotropic"    ) {
            wireMesh = new WireMesh<Symmetry::Orthotropic<>>(patternPath, 0, blendingPolySize);
        }
        else if (sym == "square"         ) {
            wireMesh = new WireMesh<Symmetry::Square<>>(patternPath, 0, blendingPolySize);
        }
        else if (sym == "triply_periodic") {
            wireMesh = new WireMesh<Symmetry::TriplyPeriodic<>>(patternPath, 0, blendingPolySize);
        }
        else if (sym == "doubly_periodic") {
            wireMesh = new WireMesh<Symmetry::DoublyPeriodic<>>(patternPath, 0, blendingPolySize);
        }
        else if (sym == "non_periodic" || sym == "3d_non_periodic") {
            wireMesh = new WireMesh<Symmetry::NonPeriodic<DEFAULT_TOL, 3>>(patternPath, 0, blendingPolySize);
        }
        else if (sym == "2d_non_periodic") {
            wireMesh = new WireMesh<Symmetry::NonPeriodic<DEFAULT_TOL, 2>>(patternPath, 0, blendingPolySize);
        }
        else throw std::runtime_error("Unknown symmetry type: " + sym);

        wireMesh->inflationGraph(params, points, edges, thicknesses, blendingParams, blendingPolyParams);

        vector<Region<Point> *> inRegions;
        vector<Region<Point> *> exceptRegions;
        extractFilteringRegions(bcondsPath, points, inRegions, exceptRegions);

        vector<Point> fixedPoints = excludedPoints(points, inRegions, exceptRegions);

        vector<int> fixedParameters;
        for (unsigned i = 0; i < fixedPoints.size(); i++) {
            vector<int> parameterIndices = wireMesh->pointToParametersIndices(fixedPoints[i], params);
            fixedParameters.insert(fixedParameters.end(), parameterIndices.begin(), parameterIndices.end());
        }

        vector<bool> solution(params.size(), false);
        for (unsigned i = 0; i < fixedParameters.size(); i++) {
            solution[fixedParameters[i]] = true;
        }

        delete wireMesh;

        return solution;
    }

    vector<bool> generateParametersMask(string patternPath, string paramsString, string bcondsPath, size_t blendingPolySize = 0, string sym = "non_periodic") {

        // Parse parameters
        auto parseParams = [](string pstring) -> vector<Real> {
            boost::trim(pstring);
            vector<string> tokens;
            boost::split(tokens, pstring, boost::is_any_of("\t "),
                         boost::token_compress_on);
            vector<Real> pvals;
            for (string &s : tokens) pvals.push_back(std::stod(s));
            return pvals;
        };

        // If requested, override the initial parameters set in the job file
        vector<double> params = parseParams(paramsString);

        return generateParametersMask(patternPath, params, bcondsPath, blendingPolySize, sym);
    }

    // Function for BoundaryPerturbationInflator
    template<size_t N>
    vector<bool> generateParametersMask(string meshPath, string bcondsPath) {
        // Original mesh
        vector<MeshIO::IOVertex>  inVertices;
        vector<MeshIO::IOElement> inElements;
        load(meshPath, inVertices, inElements, MeshIO::FMT_GUESS, MeshIO::MESH_GUESS);

        // Create inflator and initial params
        BoundaryPerturbationInflator<N> bpi(inVertices, inElements, false);
        vector<Real> params(bpi.numParameters(), 0.0);

        vector<Point> points;
        for (unsigned i=0; i<inVertices.size(); i++) {
            points.push_back(inVertices[i].point);
        }

        vector<Region<Point> *> inRegions;
        vector<Region<Point> *> exceptRegions;
        extractFilteringRegions(bcondsPath, points, inRegions, exceptRegions);

        vector<Point> fixedPoints = excludedPoints(points, inRegions, exceptRegions);

        vector<int> fixedParameters;
        for (unsigned i = 0; i < fixedPoints.size(); i++) {
            vector<int> parameterIndices = bpi.pointToParametersIndices(truncateFrom3D<VectorND<N>>(fixedPoints[i]));
            fixedParameters.insert(fixedParameters.end(), parameterIndices.begin(), parameterIndices.end());
        }

        vector<bool> solution(params.size(), false);
        for (unsigned i = 0; i < fixedParameters.size(); i++) {
            if(solution[fixedParameters[i]])
                std::cout << "Warning: parameter masked more than once!" << std::endl;
            solution[fixedParameters[i]] = true;
        }

        return solution;

    }
}

#endif //MICROSTRUCTURES_PARAMETERSMASK_H
