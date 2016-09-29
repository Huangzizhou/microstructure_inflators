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

void usage(int status, const po::options_description &visible_opts) {
    cerr << "Usage: ./isosurface_cli mesher_name pattern parameters out.msh [options]" << endl;
    cerr << "eg: ./isosurface_cli cubic pattern0746.wire \"0.25 0.5 0.25 0.25 0.5 0.25 0.25 0.25 0.25\" out.msh" << endl;
    cout << visible_opts << endl;
    exit(status);
}

po::variables_map parseCmdLine(int argc, char *argv[]) {
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("mesher", po::value<string>(), "name of mesher to use")
        ("wire", po::value<string>(), "input wire file")
        ("outMSH", po::value<string>(), "output msh file")
        ;
    po::positional_options_description p;
    p.add("mesher", 1);
    p.add("wire", 1);
    p.add("outMSH", 1);

    // Options visible in the help message.
    po::options_description visible_opts;
    visible_opts.add_options()("help,h", "Produce this help message")
                              ("disablePostprocessing,d",                    " Disable post-processing of mesher output")
                              ("dumpBaseUnitGraph,B", po::value<string>(),   " Output the base unit inflation graph to the specified path")
                              ("dumpInflationGraph,D", po::value<string>(),  " Output the inflation graph to the specified path")
                              ("dumpReplicatedGraph,R", po::value<string>(), " Output the replicated pattern graph to the specified path")
                              ("mopts,m", po::value<string>(),               " Meshing options file")
                              ("params,p", po::value<string>(),              " Pattern parameters")
                              ("nonReflectiveInflator",                      " use non-reflective inflator (reflective by default)")
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
    if (vm.count("outMSH") == 0) {
        cout << "Error: must specify all positional arguments" << endl;
        fail = true;
    }

    if (fail)
        usage(1, visible_opts);

    return vm;
}

int main(int argc, char *argv[])
{
    auto args = parseCmdLine(argc, argv);

    IsosurfaceInflator inflator(args["mesher"].as<string>(), true, args["wire"].as<string>());

    vector<Real> params(inflator.defaultParameters());
    if (args.count("params")) {
        // Split up params.
        string paramString(args["params"].as<string>());
        boost::trim(paramString);
        std::vector<string> pStrings;
        boost::split(pStrings, paramString, boost::is_any_of("\t "), boost::token_compress_on);
        params.clear();
        for (const auto &p : pStrings)
            params.push_back(std::stod(p));
    }
    else {
        cout << "Inflating default parameters: " << endl;
        for (Real p : params) cout << p << "\t";
        cout << endl;
    }

    auto &config = IsosurfaceInflatorConfig::get();
    if (args.count("dumpInflationGraph"))
        config.inflationGraphPath = args["dumpInflationGraph"].as<string>();
    if (args.count("dumpReplicatedGraph"))
        config.replicatedGraphPath = args["dumpReplicatedGraph"].as<string>();
    if (args.count("dumpBaseUnitGraph"))
        config.baseUnitGraphPath = args["dumpBaseUnitGraph"].as<string>();

    if (args.count("mopts")) inflator.meshingOptions().load(args["mopts"].as<string>());
    if (args.count("disablePostprocessing")) inflator.disablePostprocess();
    inflator.setReflectiveInflator(args.count("nonReflectiveInflator") == 0);
    inflator.inflate(params);

    MeshIO::save(args["outMSH"].as<string>(), inflator.vertices(), inflator.elements());
}
