////////////////////////////////////////////////////////////////////////////////
// Inflator_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Provides a simple command line interface to the 3D inflator.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/02/2015 02:57:00
////////////////////////////////////////////////////////////////////////////////
#include <MeshIO.hh>
#include "LinearElasticity.hh"
#include <Materials.hh>
#include <PeriodicHomogenization.hh>
#include "GlobalBenchmark.hh"

#include "Inflator.hh"

#include <vector>
#include <queue>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <random>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "PatternOptimizationJob.hh"
#include "PatternOptimizationIterate.hh"

namespace po = boost::program_options;
using namespace std;
using namespace PatternOptimization;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: Inflator_cli pattern.wire [job.opt] [options]" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("pattern", po::value<string>(), "Pattern wire file")
        ("job",  po::value<string>(), "Job file (for bounds--optional)")
        ;
    po::positional_options_description p;
    p.add("pattern", 1);
    p.add("job", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help",        "Produce this help message")
        ("cell_size,c",         po::value<double>()->default_value(5.0),           "Inflation cell size")
        ("default_thickness,t", po::value<double>()->default_value(0.5 * sqrt(2)), "default thickness (of bar diagonal)")
        ("out,o",    po::value<string>(), "Output inflated mesh")
        ("dof,d",    po::value<string>(), "dof file specifying parameters")
        ("multiDoF,D",                                    "inflate multiple DoF files (names specified on STDIN)")
        ("subdivide,S",  po::value<size_t>()->default_value(0),           "number of subdivisions to run for 3D inflator")
        ("sub_algorithm,A", po::value<string>()->default_value("simple"), "subdivision algorithm for 3D inflator (simple or loop)")
        ("max_volume,v", po::value<double>(),                             "maximum element volume parameter for wire inflator")
        ;

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    bool fail = false;
    if (vm.count("pattern") == 0) {
        cout << "Error: must specify input pattern file" << endl;
        fail = true;
    }

    set<string> subdivisionAlgorithms = {"simple", "loop"};
    if (subdivisionAlgorithms.count(vm["sub_algorithm"].as<string>()) == 0) {
        cout << "Illegal subdivision algorithm specified" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}


////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    cout << setprecision(16);
    po::variables_map args = parseCmdLine(argc, argv);

    Real defaultThickness = args["default_thickness"].as<double>();
    Real cellSize = args["cell_size"].as<double>();
    Inflator<3> inflator(args["pattern"].as<string>(),
                         cellSize, defaultThickness);
    inflator.configureSubdivision(args["sub_algorithm"].as<string>(),
                                  args["subdivide"].as<size_t>());

    if (args.count("max_volume"))
        inflator.setMaxElementVolume(args["max_volume"].as<double>());

    if (args.count("multiDoF")) {
        cout << "multiDoF mode--reading filenames from STDIN" << endl;
        string name;
        while (cin >> name) {
            inflator.inflate(name);
            cout << "inflated " << name << endl;
            assert((inflator.elements().size() > 0) &&
                   (inflator.vertices().size() > 0));
            if (args.count("output")) {
                MeshIO::save(args["out"].as<string>() + "." + name + ".msh", inflator.vertices(),
                             inflator.elements());
            }
        }
        exit(0);
    }

    // Parameter specification precedence:
    //  .dof presides over .opt, which presides over defaults.
    Real defaultOffset = 0.0;
    vector<Real> params(inflator.numParameters());
    for (size_t p = 0; p < inflator.numParameters(); ++p) {
        if (inflator.parameterType(p) == ParameterType::Thickness)
            params[p] = defaultThickness;
        else 
            params[p] = defaultOffset;
    }

    if (args.count("job")) {
        auto job = parseJobFile(args["job"].as<string>());
        params = job->initialParams;
    }

    if (args.count("dof")) {
        inflator.loadPatternDoFs(args["dof"].as<string>(), params);
    }

    // Output param types and values (if specified).
    for (size_t p = 0; p < inflator.numParameters(); ++p) {
        std::cout << "Param " << p << ": "
            << ((inflator.parameterType(p) == ParameterType::Offset) ? "offset" : "thickness")
            << ", " << params[p] << std::endl;
    }

    if (args.count("out")) {
        inflator.inflate(params);
        MeshIO::save(args["out"].as<string>(), inflator.vertices(),
                     inflator.elements());
    }
}
