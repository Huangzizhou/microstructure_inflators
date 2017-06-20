//
// Created by Davi Colli Tozoni on 6/9/17.
//
#include "IsosurfaceInflator.hh"
#include "MSHFieldWriter.hh"
#include "MeshingOptions.hh"

#include "IsosurfaceInflatorConfig.hh"

#include <iostream>
#include <boost/algorithm/string.hpp>

#include <boost/program_options.hpp>

namespace po = boost::program_options;
using namespace std;

using namespace std;

void usage(int status) {
    cerr << "Usage: ./default_parameters mesher_name pattern" << endl;
    cerr << "eg: ./isosurface_cli cubic pattern0746.wire" << endl;
    exit(status);
}

po::variables_map parseCmdLine(int argc, char *argv[]) {
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
            ("mesher", po::value<string>(), "name of mesher to use")
            ("wire", po::value<string>(), "input wire file")
            ;
    po::positional_options_description p;
    p.add("mesher", 1);
    p.add("wire", 1);

    po::options_description cli_opts;
    cli_opts.add(hidden_opts);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                options(cli_opts).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (std::exception &e) {
        cout << "Error: " << e.what() << endl << endl;
        usage(1);
    }

    return vm;
}

int main(int argc, char *argv[])
{
    auto args = parseCmdLine(argc, argv);

    IsosurfaceInflator inflator(args["mesher"].as<string>(), true, args["wire"].as<string>());

    vector<Real> params(inflator.defaultParameters());
    cout << "Inflating default parameters: " << endl;
    for (Real p : params) cout << p << "\t";
    cout << endl;
}

