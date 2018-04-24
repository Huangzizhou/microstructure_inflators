////////////////////////////////////////////////////////////////////////////////
// ParamRegion_cli.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Attempt to visualize the valid parameter subspace.
*/ 
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  01/05/2015 04:49:48
////////////////////////////////////////////////////////////////////////////////
#include <MeshFEM/MeshIO.hh>
#include "LinearElasticity.hh"
#include <MeshFEM/Materials.hh>
#include <MeshFEM/PeriodicHomogenization.hh>
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
    cout << "Usage: ParamRegion_cli pattern.wire outPrefix p0 p1 [options]" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("pattern", po::value<string>(), "Pattern wire file")
        ("out",  po::value<string>(), "PGM output filename prefix")
        ("p0",   po::value<int>(), "Parameter 0 (vertical axis)")
        ("p1",   po::value<int>(), "Parameter 1 (horizontal axis)")
        ;
    po::positional_options_description p;
    p.add("pattern", 1);
    p.add("out", 1);
    p.add("p0", 1);
    p.add("p1", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help",        "Produce this help message")
        ("sampleSize,s",           po::value<int>()->default_value(100),              "samples along each parameter")
        ("cell_size,c",            po::value<double>()->default_value(5.0),           "Inflation cell size")
        ("isotropicParameters,I",  "Use isotropic DoFs")
        ("vertexThickness,V",      "Use vertex thickness instead of edge thickness")
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
    if (vm.count("p1") == 0) {
        cout << "Error: must specify input pattern, output image prefix, and two parameter indices" << endl;
        fail = true;
    }

    if (vm.count("pattern") == 0) {
        cout << "Error: must specify input pattern file" << endl;
        fail = true;
    }

    if (vm["sampleSize"].as<int>() < 2) {
        cout << "Must use at least two samples per param" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

struct Bound {
    Bound() { set(0, 0, 0);  }
    void set(Real l, Real u, Real d) { lower = l; upper = u; defaultVal = d; }
    Real midpoint() const { return 0.5 * (lower + upper); } 
    Real width() const { return upper - lower; }
    Real lower, upper, defaultVal;
};

void setDefault(vector<Real> &params, const vector<Bound> &bounds) {
    for (size_t i = 0; i < params.size(); ++i)
        params[i] = bounds.at(i).defaultVal;
}

void print(ostream &os, const vector<Real> &params) {
    for (Real p : params) {
        os << p << "\t";
    }
    os << endl;
}

////////////////////////////////////////////////////////////////////////////////
/*! Program entry point
//  @param[in]  argc    Number of arguments
//  @param[in]  argv    Argument strings
//  @return     status  (0 on success)
*///////////////////////////////////////////////////////////////////////////////
int main(int argc, const char *argv[])
{
    cout << setprecision(20);
    po::variables_map args = parseCmdLine(argc, argv);

    Real cellSize = args["cell_size"].as<double>();
    bool isotropicParameters = args.count("isotropicParameters");
    bool vertexThickness = args.count("vertexThickness");
    Inflator<3> inflator(args["pattern"].as<string>(), cellSize,
                         0.5 * sqrt(2), isotropicParameters,
                         vertexThickness);

    size_t nParams = inflator.numParameters();
    vector<Bound> bounds(nParams);
    vector<Real> params(nParams);

    // Parameter specification precedence:
    //  .dof presides over .opt, which presides over defaults.
    for (size_t p = 0; p < inflator.numParameters(); ++p) {
        if (inflator.parameterType(p) == ParameterType::Thickness)
            bounds[p].set(0.4, 1.0, 0.4 * sqrt(2));
        else 
            bounds[p].set(-0.40, 0.40, 0.0);
    }
    
    size_t nSamples = args["sampleSize"].as<int>();
    vector<unsigned char> imgBuffer(nSamples * nSamples);

    size_t p0 = args["p0"].as<int>();
    size_t p1 = args["p1"].as<int>();

    if ((p0 >= nParams) || (p1 >= nParams))
        throw runtime_error("Invalid slice parameter indices");

    string outPrefix = args["out"].as<string>();
    string outName = outPrefix + "." + to_string(p0) + "." + to_string(p1);
    // ofstream outFile(outName + ".txt");

    setDefault(params, bounds);
    for (size_t v0 = 0; v0 < nSamples; ++v0) {
        params[p0] = bounds[p0].lower + (Real(v0) * bounds[p1].width()) / (nSamples - 1);
        for (size_t v1 = 0; v1 < nSamples; ++v1) {
            params[p1] = bounds[p1].lower + (Real(v1) * bounds[p1].width()) / (nSamples - 1);
            print(cout, params);
            unsigned char val = 255;
            try { inflator.inflate(params); }
            catch (std::exception &e) { val = 0; cerr << e.what() << endl;}
            imgBuffer[v0 * nSamples + v1] = val;
        }
    }

    // Output PGM
    ofstream pgmFile(outName + ".pgm");
    if (!pgmFile.is_open()) throw runtime_error("Couldn't open output image");
    pgmFile << "P5" << endl << nSamples << "\t" << nSamples << endl << 255 << endl;
    pgmFile.write((char *) &imgBuffer[0], imgBuffer.size());
    pgmFile.close();
}
