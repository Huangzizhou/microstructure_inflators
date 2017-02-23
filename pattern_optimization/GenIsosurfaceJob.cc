////////////////////////////////////////////////////////////////////////////////
// GenIsosurfaceJob.cc
////////////////////////////////////////////////////////////////////////////////
/*! @file
//      Create a pattern optimization job intended to be used with the
//      isosurface inflator. This is particularly useful for generating sane
//      initial parameters (from the defaults) and from converting from offset
//      bounds to translation bounds:
//
//          The new isosurface inflator uses translations as variables because these
//          allow the printability constraints to be formulated more transparently.
//
//          However, often it is desirable to bound the offset of vertices from
//          known good positions (e.g. positions initialized from topology
//          enumeration). This utility will express such an offset bound as a
//          separate bound on each translation variable.
*/
//  Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
//  Company:  New York University
//  Created:  07/04/2016 17:47:39
////////////////////////////////////////////////////////////////////////////////
#include <boost/program_options.hpp>

#include <iostream>
#include <cstdlib>
#include <string>
#include <cassert>
#include <stdexcept>

#include <Future.hh>
#include "PatternOptimizationJob.hh"
#include <cmath>

#include <boost/algorithm/string.hpp>

#include "../isosurface_inflator/WireMesh.hh"

namespace po = boost::program_options;
using namespace std;

void usage(int exitVal, const po::options_description &visible_opts) {
    cout << "Usage: GenIsosurfaceJob [options] topology.wire" << endl;
    cout << visible_opts << endl;
    exit(exitVal);
}

po::variables_map parseCmdLine(int argc, const char *argv[])
{
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("topology", po::value<string>(), "Topology (line mesh)")
        ;
    po::positional_options_description p;
    p.add("topology", 1);

    po::options_description visible_opts;
    visible_opts.add_options()("help",        "Produce this help message")
        ("offsetBounds,o",      po::value<string>(),                             "offset bounds specifier (lower,upper)")
        ("translationBounds,t", po::value<string>(),                             "translation bounds specifier (lower,upper)")
        ("radiusBounds,r",      po::value<string>()->default_value("0.04,0.2"),  "radius bounds specifier (lower,upper)")
        ("blendingBounds,b",    po::value<string>()->default_value("0.005,0.2"), "blending bounds specifier (lower,upper)")
        ("elasticityTensor,e",  po::value<string>()->default_value("1,0"),       "target tensor specifier (Young,Poisson)")
        ("initialParams,p",     po::value<string>(),                             "initial parameters (optional)")
        ("parameterConstraints,c", po::value<string>(),                          "parameter constraint expressions (semicolon-separated, optional)")
        ;

    po::options_description cli_opts;
    cli_opts.add(visible_opts).add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                  options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    bool fail = false;

    if (vm.count("topology") == 0) {
        cout << "Error: must specify pattern topology" << endl;
        fail = true;
    }

    if (vm.count("translationBounds") && vm.count("offsetBounds")) {
        cout << "Error: must not specify both translationBounds and offsetBounds" << endl;
        fail = true;
    }

    if (fail || vm.count("help"))
        usage(fail, visible_opts);

    return vm;
}

vector<Real> parseVecArg(const string &arg, size_t N) {
    vector<string> components;
    boost::split(components, arg, boost::is_any_of(","));
    if (components.size() != N) throw runtime_error("Invalid number of components in argument " + arg);

    vector<Real> result;
    for (const string &c : components)
        result.push_back(stod(c));
    return result;
}

int main(int argc, const char *argv[])
{
    auto args = parseCmdLine(argc, argv);

    // Infer dimension from pattern
    vector<MeshIO::IOVertex> vertices;
    vector<MeshIO::IOElement> elements;
    MeshIO::load(args.at("topology").as<string>(), vertices, elements);

    Real zmag = 0;
    for (const auto &v : vertices)
        zmag = max(zmag, std::abs(v[2]));
    size_t dim = (zmag < 1e-2) ? 2 : 3;

    WireMesh<ThicknessType::Vertex, Symmetry::Orthotropic<>> wm(vertices, elements);

    vector<Real> targetModuli = parseVecArg(args["elasticityTensor"].as<string>(), 2);

    unique_ptr<PatternOptimization::JobBase> job;
    if (dim == 2) {
        auto job2D = Future::make_unique<PatternOptimization::Job<2>>();
        job2D->targetMaterial.setIsotropic(targetModuli[0], targetModuli[1]);
        job = move(job2D);
    }
    else {
        auto job3D = Future::make_unique<PatternOptimization::Job<3>>();
        job3D->targetMaterial.setIsotropic(targetModuli[0], targetModuli[1]);
        job = move(job3D);
    }

    if (args.count("offsetBounds")) {
        vector<Real> offsetBds = parseVecArg(args["offsetBounds"].as<string>(), 2);
        auto defaultPositions = wm.defaultPositionParams();

        for (size_t p = 0; p < defaultPositions.size(); ++p) {
            // Position parameters should be first in the isosurface inflator
            assert(wm.isPositionParam(p));
            job->varLowerBounds.emplace(p, defaultPositions[p] + offsetBds[0]);
            job->varUpperBounds.emplace(p, defaultPositions[p] + offsetBds[1]);
        }
    }

    // Set translation bounds, which will be ignored by the optimizer if
    // offsetBounds introduced per-variable bounds above
    job->translationBounds = {0.1, 0.8};
    if (args.count("translationBounds"))
        job->translationBounds = parseVecArg(args["translationBounds"].as<string>(), 2);

    job->radiusBounds   = parseVecArg(args[  "radiusBounds"].as<string>(), 2);
    job->blendingBounds = parseVecArg(args["blendingBounds"].as<string>(), 2);

    if (args.count("parameterConstraints")) {
        boost::split(job->parameterConstraints,
                args["parameterConstraints"].as<string>(), boost::is_any_of(";"));
    }

    job->initialParams = wm.defaultParameters();
    if (args.count("initialParams")) {
        job->initialParams = parseVecArg(args["initialParams"].as<string>(),
                                         job->initialParams.size());
    }

    job->writeJobFile(cout);

    return 0;
}
