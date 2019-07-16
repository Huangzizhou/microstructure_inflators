////////////////////////////////////////////////////////////////////////////////
// ParametersMask.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Given a pattern, a set of parameters and a set of boundary conditions,
//    decide which variables should not change during optimization.
//      This is done by returning vector of booleans with 'false' for variables
//    that can move and 'true' for variables that should be filtered out of the
//    optimization.
//
*/
//  Author:  Davi Colli Tozoni (dctozoni) davi.tozoni@nyu.edu
//  Company:  New York University
//  Created:  2/13/18
////////////////////////////////////////////////////////////////////////////////

#ifndef PARAMETERSMASK_HH
#define PARAMETERSMASK_HH

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
    // Function that loads json file, creating regions that should be included or excluded from optimization
     void jsonToRegions(std::string jsonPath, const BBox<Point> &bbox, vector<Region<Point> *> &inRegions, vector<Region<Point> *> &exceptRegions, vector<Region<Point> *> &glueRegions) {

        // reading JSON file
        ifstream input(jsonPath);
        json j;
        input >> j;

        json regions = j["regions"];
        // iterate the array
        for (json::iterator it = regions.begin(); it != regions.end(); ++it) {
            Region<Point> * newRegion;

            Point minCorner;
            Point maxCorner;

            json current_region = *it;

            json box = current_region["box"];
            json box_percentage = current_region["box%"];
            json path = current_region["path"];

            if (box.size() > 0) {
                vector<float> min_corner = box["minCorner"];
                vector<float> max_corner = box["maxCorner"];

                minCorner << 0.0, 0.0, 0.0;
                maxCorner << 0.0, 0.0, 0.0;

                for (size_t i=0; i<min_corner.size(); i++) {
                    minCorner[i] = min_corner[i];
                }

                for (size_t i=0; i<max_corner.size(); i++) {
                    maxCorner[i] = max_corner[i];
                }

                newRegion = new BBox<Point>(minCorner, maxCorner);
            }
            if (box_percentage.size() > 0) {
                // Find the bounding box
                Eigen::Vector3d m = bbox.minCorner;
                Eigen::Vector3d M = bbox.maxCorner;

                vector<float> min_corner = box_percentage["minCorner"];
                vector<float> max_corner = box_percentage["maxCorner"];

                Eigen::Array3d min_corner_array;
                min_corner_array << 0.0, 0.0, 0.0;
                Eigen::Array3d max_corner_array;
                max_corner_array << 0.0, 0.0, 0.0;

                for (size_t i=0; i<min_corner.size(); i++) {
                    min_corner_array[i] = min_corner[i];
                }

                for (size_t i=0; i<max_corner.size(); i++) {
                    max_corner_array[i] = max_corner[i];
                }

                // scale
                double w = M(0) - m(0);
                double h = M(1) - m(1);
                double d = M(2) - m(2);
                min_corner_array(0) *= w;
                min_corner_array(1) *= h;
                min_corner_array(2) *= d;
                max_corner_array(0) *= w;
                max_corner_array(1) *= h;
                max_corner_array(2) *= d;

                // translate to min corner
                min_corner_array += m.array();
                max_corner_array += m.array();

                minCorner = min_corner_array.matrix();
                maxCorner = max_corner_array.matrix();

                newRegion = new BBox<Point>(minCorner, maxCorner);
            }
            if (path.size() > 0) {
                // Always 2D
                std::vector<Point> points;
                for (std::vector<double> point : path) {
                    Point eigenPoint;
                    for (size_t i=0; i<point.size(); i++) {
                        eigenPoint[i] = point[i];
                    }
                    eigenPoint[2] = 0.0;
                    points.push_back(eigenPoint);
                }

                newRegion = new PathRegion<Point>(points);
            }

            string typeString = current_region["type"];
            if (typeString.compare("dirichlet") == 0 || typeString.compare("force") == 0 || typeString.compare("zero") == 0 || typeString.compare("contact") == 0) {
                exceptRegions.push_back(newRegion);
            }
            else if (typeString.compare("optimization") == 0) {
                inRegions.push_back(newRegion);
            }
            else if (typeString.compare("glue") == 0 || typeString.compare("fracture") == 0 ) {
                glueRegions.push_back(newRegion);
            }
            else {
                std::cerr << "Region of type " + typeString + " not expected." << std::endl;
                throw std::runtime_error("Region of type " + typeString + " not expected.");
            }
        }
    }

    // Decide which points should be excluded from optimization
    vector<Point> excludedPoints(vector<Point> points, vector<Region<Point> *> &inRegions, vector<Region<Point> *> &exceptRegions) {
        set<int> in;

        // saving points that are inside 'in' regions
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
            // In case there are no 'in' regions, consider whole region as in.
            for (size_t i = 0; i < points.size(); i++) {
                in.insert(i);
            }
        }

        // remove points from 'in' regions that are also inside except regions
        for (auto region : exceptRegions) {
            for (size_t i = 0; i < points.size(); i++) {
                if (region->containsPoint(points[i])) {
                    in.erase(i);
                }
            }
        }

        // create resulting vector which is complement of in set
        vector<Point> result;
        for (size_t i = 0; i < points.size(); i++) {
            if (!in.count(i))
                result.push_back(points[i]);
        }

        return result;
    }

    // Extract filtering regions from file
    void extractFilteringRegions(string bcondsPath, vector<Point> points, vector<Region<Point> *> &inRegions, vector<Region<Point> *> &exceptRegions, vector<Region<Point> *> &glueRegions) {
        if (!bcondsPath.empty()) {
            BBox<Point> bb(points);

            jsonToRegions(bcondsPath, bb, inRegions, exceptRegions, glueRegions);
        }
    }

    // Generate parameters mask when using truss like structures
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
        vector<Region<Point> *> glueRegions;
        extractFilteringRegions(bcondsPath, points, inRegions, exceptRegions, glueRegions);

        // find which points should be excluded
        vector<Point> fixedPoints = excludedPoints(points, inRegions, exceptRegions);

        // transform points into parameters
        vector<int> fixedParameters;
        for (unsigned i = 0; i < fixedPoints.size(); i++) {
            vector<int> parameterIndices = wireMesh->pointToParametersIndices(fixedPoints[i], params);
            fixedParameters.insert(fixedParameters.end(), parameterIndices.begin(), parameterIndices.end());
        }

        // create filtering vector
        vector<bool> solution(params.size(), false);
        for (unsigned i = 0; i < fixedParameters.size(); i++) {
            solution[fixedParameters[i]] = true;
        }

        delete wireMesh;

        return solution;
    }

    // Generate parameters mask when using truss like structures (when parameters are passed as string)
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

    // Generate parameters mask when using boundary perturbation inflator (meaning that each coordinate of each boundary
    // vertex could be a variable in the optimization)
    template<size_t N>
    vector<bool> generateParametersMask(const vector<MeshIO::IOVertex> &inVertices, const vector<MeshIO::IOElement> &inElements, string bcondsPath) {
        // Create inflator and initial params
        BoundaryPerturbationInflator<N> bpi(inVertices, inElements, false);
        vector<Real> params(bpi.numParameters(), 0.0);

        vector<Point> points;
        for (unsigned i=0; i<inVertices.size(); i++) {
            points.push_back(inVertices[i].point);
        }

        vector<Region<Point> *> inRegions;
        vector<Region<Point> *> exceptRegions;
        vector<Region<Point> *> glueRegions;
        extractFilteringRegions(bcondsPath, points, inRegions, exceptRegions, glueRegions);

        // find which points should be excluded
        vector<Point> fixedPoints = excludedPoints(points, inRegions, exceptRegions);

        // transform points into parameters
        vector<int> fixedParameters;
        for (unsigned i = 0; i < fixedPoints.size(); i++) {
            vector<int> parameterIndices = bpi.pointToParametersIndices(truncateFrom3D<VectorND<N>>(fixedPoints[i]));
            fixedParameters.insert(fixedParameters.end(), parameterIndices.begin(), parameterIndices.end());
        }

        vector<bool> solution(params.size(), false);
        for (unsigned i = 0; i < fixedParameters.size(); i++) {
            if(solution[fixedParameters[i]]) {
                std::cout << "Warning: parameter masked more than once!" << std::endl;
                //int param = fixedParameters[i];
                //std::cout << "Param: " << param << std::endl;
                //std::cout << "Point: " << bpi.parameterIndexToPoint(param) << std::endl;
            }

            solution[fixedParameters[i]] = true;
        }

        return solution;
    }

    // Generate parameters mask when using boundary perturbation inflator (meaning that each coordinate of each boundary
    // vertex could be a variable in the optimization)
    template<size_t N>
    vector<bool> generateParametersMask(string meshPath, string bcondsPath) {
        // Original mesh
        vector<MeshIO::IOVertex>  inVertices;
        vector<MeshIO::IOElement> inElements;
        load(meshPath, inVertices, inElements, MeshIO::FMT_GUESS, MeshIO::MESH_GUESS);

        return generateParametersMask<N>(inVertices, inElements, bcondsPath);
    }


    // Generate map of parameters that are dependent on others. In the optimization, the corresponding vertices
    // will always move together. Mostly used for contact optimization
    template<size_t N>
    vector<int> generateGluedParametersMap(const vector<MeshIO::IOVertex> &inVertices, const vector<MeshIO::IOElement> &inElements, string bcondsPath) {
        // Create inflator and initial params
        BoundaryPerturbationInflator<N> bpi(inVertices, inElements, false);
        vector<Real> params(bpi.numParameters(), 0.0);
        vector<int> parameterToGluedParameter(bpi.numParameters(), -1);

        vector<Point> points;
        for (unsigned i=0; i<inVertices.size(); i++) {
            points.push_back(inVertices[i].point);
        }

        vector<Region<Point> *> inRegions;
        vector<Region<Point> *> exceptRegions;
        vector<Region<Point> *> glueRegions;
        extractFilteringRegions(bcondsPath, points, inRegions, exceptRegions, glueRegions);

        if (glueRegions.size() > 0) {
            // Construct glued params vector
            vector<Point> gluedPoints = excludedPoints(points, inRegions, glueRegions);
            parameterToGluedParameter.assign(params.size(), -1);
            for (unsigned i = 0; i < gluedPoints.size(); i++) {
                vector<int> parameterIndices = bpi.pointToParametersIndices(truncateFrom3D<VectorND<N>>(gluedPoints[i]));

                if (parameterIndices.size() == 2 * N) {
                    for (size_t d = 0; d < N; d++) {
                        parameterToGluedParameter[parameterIndices[d]] = parameterIndices[N + d];
                        parameterToGluedParameter[parameterIndices[N + d]] = parameterIndices[d];
                    }

                }
            }
        }

        return parameterToGluedParameter;
    }

    // Generate map of parameters that are dependent on others. In the optimization, the corresponding vertices
    // will always move together. Mostly used for contact optimization
    template<size_t N>
    vector<int> generateGluedParametersMap(string meshPath, string bcondsPath) {
        // Original mesh
        vector<MeshIO::IOVertex>  inVertices;
        vector<MeshIO::IOElement> inElements;
        load(meshPath, inVertices, inElements, MeshIO::FMT_GUESS, MeshIO::MESH_GUESS);

        return generateGluedParametersMap<N>(inVertices, inElements, bcondsPath);
    }
}

#endif //PARAMETERSMASK_HH
