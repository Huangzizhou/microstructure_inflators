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
#include <MeshFEM/MeshIO.hh>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include "ParametersMask.hh"

using namespace std;

namespace po = boost::program_options;

using Point = Point3<double>;


void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: GenerateParametersMask --params '0.0 1.0 ...' --pattern octa_cell.obj --boundaryConditions bc.json" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}


po::variables_map parseCmdLine(int argc, const char *argv[]) {
    po::options_description patternOptions;
    patternOptions.add_options()
            ("pattern,p", po::value<string>(),
             "Pattern wire mesh (.obj|wire), or initial mesh for BoundaryPerturbationInflator")
            ("params", po::value<string>(), "Initial params (overrides those specified in job file).");

    po::options_description simulationOptions;
    simulationOptions.add_options()
            ("boundaryConditions,b", po::value<string>(), "boundary conditions");

    po::options_description generalOptions;
    generalOptions.add_options()
            ("help,h", "Produce this help message");

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


void execute(po::variables_map args) {
    cout << "Executing..." << endl;

    vector<bool> solution = ParametersMask::generateParametersMask(args["pattern"].as<string>(), args["params"].as<string>(), args["boundaryConditions"].as<string>());

    // print solution
    for (unsigned i = 0; i < solution.size(); i++) {
        if (solution[i])
            cout << 1 << " ";
        else
            cout << 0 << " ";
    }
    cout << endl;
}


int main(int argc, const char *argv[]) {
    po::variables_map args = parseCmdLine(argc, argv);

    execute(args);

    return 0;
}


