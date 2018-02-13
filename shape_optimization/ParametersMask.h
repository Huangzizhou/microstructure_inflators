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
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <LinearElasticity.hh>
#include <BoundaryConditions.hh>
#include <Future.hh>
#include <WireMesh.hh>
#include <Symmetry.hh>

using namespace std;

namespace po = boost::program_options;

using Point = Point3<double>;

namespace ParametersMask {


    vector<Point> pointsInsideBCs(vector<Point> points, vector<BBox<Point>> regions) {
        vector<Point> result;

        // marking which points are filtered and which are not
        for (auto region : regions) {
            for (int i = 0; i < points.size(); i++) {
                if (region.containsPoint(points[i])) {
                    result.push_back(points[i]);
                } else {
                    // so far, nothing
                }
            }
        }

        return result;
    }

    vector<BBox<Point>> extractFilteringRegions(string bcondsPath, vector<Point> points) {
        vector<BBox<Point>> filteringRegions;
        if (!bcondsPath.empty()) {
            BBox<Point> bb(points);
            auto dim = bb.dimensions();
            if ((std::abs(dim[0]) < 1e-6) || (std::abs(dim[1]) < 1e-6))
                throw std::runtime_error("Degenerate pattern");
            const int N = std::abs(dim[2]) <= 1e-6 ? 2 : 3;

            if (N == 2) {
                bool no_rigid_motion;
                Point2D bb_min(bb.minCorner(0), bb.minCorner(1));
                Point2D bb_max(bb.maxCorner(0), bb.maxCorner(1));
                BBox<Point2D> bb_2D(bb_min, bb_max);

                std::vector<CondPtr<2>> boundary_conditions = readBoundaryConditions<2>(bcondsPath, bb_2D,
                                                                                        no_rigid_motion);
                for (auto boundary : boundary_conditions) {
                    BBox<Point2D> new_region_2D = boundary.get()->region;
                    Point min_corner(new_region_2D.minCorner(0), new_region_2D.minCorner(1), 0.0);
                    Point max_corner(new_region_2D.maxCorner(0), new_region_2D.maxCorner(1), 0.0);;

                    BBox<Point> new_region(min_corner, max_corner);
                    filteringRegions.push_back(new_region);
                }
            } else {
                bool no_rigid_motion;
                std::vector<CondPtr<3>> boundary_conditions = readBoundaryConditions<3>(bcondsPath, bb,
                                                                                        no_rigid_motion);
                for (auto boundary : boundary_conditions) {
                    BBox<Point> new_region = boundary.get()->region;
                    filteringRegions.push_back(new_region);
                }
            }
        }

        return filteringRegions;
    }

    vector<bool> generateParametersMask(string patternPath, string paramsString, string bcondsPath) {
        WireMesh<Symmetry::NonPeriodic<DEFAULT_TOL>> wireMesh(patternPath);

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

        vector<Point3<double>> points;
        vector<pair<size_t, size_t>> edges;
        vector<double> thicknesses;
        vector<double> blendingParams;

        wireMesh.inflationGraph(params, points, edges, thicknesses, blendingParams);

        vector<BBox<Point>> filteringRegions = extractFilteringRegions(bcondsPath, points);

        vector<Point> fixedPoints = pointsInsideBCs(points, filteringRegions);

        vector<int> fixedParameters;
        for (unsigned i = 0; i < fixedPoints.size(); i++) {
            vector<int> parameterIndices = wireMesh.pointToParametersIndices(fixedPoints[i]);
            fixedParameters.insert(fixedParameters.end(), parameterIndices.begin(), parameterIndices.end());
        }

        vector<bool> solution(params.size(), false);
        for (unsigned i = 0; i < fixedParameters.size(); i++) {
            solution[fixedParameters[i]] = true;
        }

        return solution;
    }
}

#endif //MICROSTRUCTURES_PARAMETERSMASK_H
