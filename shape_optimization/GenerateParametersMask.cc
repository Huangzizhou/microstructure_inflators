////////////////////////////////////////////////////////////////////////////////
// GenerateParametersMask.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Given a pattern, a set of parameters and a set of boundary conditions,
//    decide which variables should not change during optimization.
//
*/
//  Author:  Davi Colli Tozoni (dctozoni) davi.tozoni@nyu.edu
//  Company:  New York University
//  Created:  2/12/18
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include <LinearElasticity.hh>
#include <BoundaryConditions.hh>
#include <Future.hh>
#include <vector>
#include <iostream>
#include <WireMesh.hh>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <Symmetry.hh>

using namespace std;

namespace po = boost::program_options;

using Point = Point3<double>;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: GenerateParametersMask --params '0.0 1.0 ...' --pattern octa_cell.obj --boundaryConditions bc.json" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description patternOptions;
    patternOptions.add_options()
            ("pattern,p",    po::value<string>(),                              "Pattern wire mesh (.obj|wire), or initial mesh for BoundaryPerturbationInflator")
            ("params",       po::value<string>(),                              "Initial params (overrides those specified in job file).")
            ;

    po::options_description simulationOptions;
    simulationOptions.add_options()
            ("boundaryConditions,b", po::value<string>(),                    "boundary conditions")
            ;

    po::options_description generalOptions;
    generalOptions.add_options()
            ("help,h",                                    "Produce this help message")
            ;

    po::options_description visibleOptions;
    visibleOptions.add(patternOptions).add(simulationOptions).add(generalOptions);

    po::options_description cli_opts;
    cli_opts.add(visibleOptions);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                options(cli_opts).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visibleOptions);
    }

    bool fail = false;
    if (vm.count("pattern") == 0 || vm.count("params") == 0 || vm.count("boundaryConditions") == 0) {
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visibleOptions);

    return vm;
}

vector<Point> pointsInsideBCs(vector<Point> points, vector<BBox<Point>> regions) {
    vector<Point> result;

    // marking which points are filtered and which are not
    for (auto region : regions) {
        for (int i=0; i<points.size(); i++) {
            if (region.containsPoint(points[i])) {
                result.push_back(points[i]);
            }
            else {
                // so far, nothing
            }
        }
    }

    return result;
}

void execute(po::variables_map args) {
    cout << "Executing..." << endl;

    WireMesh<Symmetry::NonPeriodic<DEFAULT_TOL>> wireMesh(args["pattern"].as<string>());

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
    vector<double> params = parseParams(args["params"].as<string>());

    vector<Point3<double>> points;
    vector<pair<size_t, size_t>> edges;
    vector<double> thicknesses;
    vector<double> blendingParams;

    wireMesh.inflationGraph(params, points, edges, thicknesses, blendingParams);

    // Now, read regions
    vector<BBox<Point>> filteringRegions;
    string bconds_path = args["boundaryConditions"].as<string>();
    if (!bconds_path.empty()) {
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

            std::vector<CondPtr<2>> boundary_conditions = readBoundaryConditions<2>(bconds_path, bb_2D, no_rigid_motion);
            for (auto boundary : boundary_conditions) {
                BBox<Point2D> new_region_2D = boundary.get()->region;
                Point min_corner(new_region_2D.minCorner(0), new_region_2D.minCorner(1), 0.0);
                Point max_corner(new_region_2D.maxCorner(0), new_region_2D.maxCorner(1), 0.0);;

                BBox<Point> new_region(min_corner, max_corner);
                filteringRegions.push_back(new_region);
            }
        }
        else {
            bool no_rigid_motion;
            std::vector<CondPtr<3>> boundary_conditions = readBoundaryConditions<3>(bconds_path, bb, no_rigid_motion);
            for (auto boundary : boundary_conditions) {
                BBox<Point> new_region = boundary.get()->region;
                filteringRegions.push_back(new_region);
            }
        }
    }

    vector<Point> fixedPoints = pointsInsideBCs(points, filteringRegions);

    vector<int> fixedParameters;
    for (unsigned i = 0; i < fixedPoints.size(); i++) {
        cout << "Filtered out: " << fixedPoints[i] << endl;

        vector<int> parameterIndices = wireMesh.pointToParametersIndices(fixedPoints[i]);
        fixedParameters.insert(fixedParameters.end(), parameterIndices.begin(), parameterIndices.end());
    }

    vector<int> solution(params.size(), 0);
    for (unsigned i = 0; i < fixedParameters.size(); i++) {
        solution[fixedParameters[i]] = 1;
    }


    // print solution
    for (unsigned i = 0; i < solution.size(); i++) {
        cout << solution[i] << " ";
    }
    cout << endl;
}


int main(int argc, const char *argv[]) {
    po::variables_map args = parseCmdLine(argc, argv);

    execute(args);

    return 0;
}


