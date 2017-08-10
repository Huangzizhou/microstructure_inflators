#include "IsosurfaceInflator.hh"
#include <boost/algorithm/string.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;

void usage(int status) {
    cerr << "Usage: ./default_parameters mesher_name pattern" << endl;
    cerr << "eg: ./default_parameters cubic pattern0746.wire" << endl;
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

    // Options visible in the help message.
    po::options_description visible_opts;
    visible_opts.add_options()("help,h", "Produce this help message")
            ("inflation_graph_radius", po::value<size_t>()->default_value(2),   "Number of edges to traverse outward from the symmetry cell when building the inflation graph")
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
        usage(1);
    }

    bool fail = false;
    if (vm.count("wire") == 0) {
        cout << "Error: must specify mesher and pattern" << endl;
        fail = true;
    }

    if (fail)
        usage(1);

    return vm;
}

int main(int argc, char *argv[])
{
    auto args = parseCmdLine(argc, argv);

    IsosurfaceInflator inflator(args["mesher"].as<string>(), true, args["wire"].as<string>(),
                                args["inflation_graph_radius"].as<size_t>());

    vector<Real> params(inflator.defaultParameters());

    cout << "Inflating Default parameters: " << endl;
    for (Real p : params) cout << p << "\t";
    cout << endl;

}
