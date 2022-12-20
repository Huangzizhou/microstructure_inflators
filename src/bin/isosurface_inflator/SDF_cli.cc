#include <isosurface_inflator/WireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/IsosurfaceInflatorConfig.hh>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <openvdb/openvdb.h>

#include <iostream>
#include <functional>

namespace po = boost::program_options;
using namespace std;

template<class WMesh>
void execute(const string &wirePath, string paramString, const string &igraphPath) {
    WMesh wmesh(wirePath);
    PatternSignedDistance<Real, WMesh> psd(wmesh);

    vector<Real> params(wmesh.defaultParameters());
    if (paramString != "") {
        // Split up params.
        boost::trim(paramString);
        std::vector<string> pStrings;
        boost::split(pStrings, paramString, boost::is_any_of("\t "), boost::token_compress_on);
        params.clear();
        for (const auto &p : pStrings)
            params.push_back(std::stod(p));
    }
    else {
        cout << "Using default parameters: " << endl;
        for (Real p : params) cout << p << "\t";
        cout << endl;
    }

    if (igraphPath != "") wmesh.saveInflationGraph(igraphPath, params);

    psd.setParameters(params, Eigen::Matrix3d::Identity());

    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create();
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    
    const auto bbox = psd.boundingBox();
    Point3d minCorner = bbox.minCorner;
    Point3d maxCorner = bbox.maxCorner;
    Eigen::Vector3i res(100, 100, 100);
    for (int i = 0; i <= res(0); i++)
    {
        const double x = minCorner(0) + (maxCorner(0) - minCorner(0)) * ((double)i / res(0)); 
        for (int j = 0; j <= res(1); j++)
        {
            const double y = minCorner(1) + (maxCorner(1) - minCorner(1)) * ((double)j / res(1)); 
            for (int k = 0; k <= res(2); k++)
            {
                const double z = minCorner(2) + (maxCorner(2) - minCorner(2)) * ((double)k / res(2)); 

                Point3d evalPt(x, y, z);
                openvdb::Coord xyz(i, j, k);
                accessor.setValue(xyz, psd.signedDistance(evalPt));
            }
        }
    }

    grid->setName("level_set");
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    openvdb::io::File file("implicit.vdb");

    openvdb::GridPtrVec(grids);
    grids.push_back(grid);
    file.write(grids);
    file.close();
}

map<string, function<void(const string &, string, const string &)>> execImpls =
    {{          "cubic", execute<WireMesh<Symmetry::Cubic<>>         >},
     {    "orthotropic", execute<WireMesh<Symmetry::Orthotropic<>>   >},
     {"triply_periodic", execute<WireMesh<Symmetry::TriplyPeriodic<>>>}};

void usage(int status, const po::options_description &visible_opts) {
    cerr << "Usage: ./SDF_cli pattern_symmetry wire.obj -p parameters [options]" << endl;
    cerr << "eg: ./SDF_cli cubic pattern0746.wire -p \"0.25 0.5 0.25 0.25 0.5 0.25 0.25 0.25 0.25\"" << endl;
    cout << visible_opts << endl;
    exit(status);
}

po::variables_map parseCmdLine(int argc, char *argv[]) {
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("symmetry",  po::value<string>(), "Pattern symmetry (triply_periodic,  orthotropic,  cubic)")
        ("wire"    ,  po::value<string>(), "input wire file")
        ;
    po::positional_options_description p;
    p.add("symmetry", 1);
    p.add("wire", 1);

    // Options visible in the help message.
    po::options_description visible_opts;
    visible_opts.add_options()("help,h", "Produce this help message")
                              ("dumpInflationGraph,D", po::value<string>()->default_value(""),  " Output the inflation graph to the specified path")
                              ("params,p", po::value<string>(),              " Pattern parameters")
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
        cerr << "Error: " << e.what() << endl << endl;
        usage(1, visible_opts);
    }

    bool fail = false;
    {
        auto symmetry = vm["symmetry"].as<string>();
        if (execImpls.count(symmetry) == 0) {
            cerr << "Error: invalid symmetry option '" << symmetry << "'" << endl;
            fail = true;
        }
    }

    if (fail)
        usage(1, visible_opts);

    return vm;
}

int main(int argc, char *argv[])
{
    auto args = parseCmdLine(argc, argv);

    auto &config = IsosurfaceInflatorConfig::get();
    if (args.count("dumpInflationGraph"))
        config.inflationGraphPath = args["dumpInflationGraph"].as<string>();
    if (args.count("dumpReplicatedGraph"))
        config.replicatedGraphPath = args["dumpReplicatedGraph"].as<string>();

    auto &exec = execImpls.at(args["symmetry"].as<string>());

    string paramString;
    if (args.count("params")) paramString = args["params"].as<string>();

    exec(args["wire"].as<string>(), paramString,
         args["dumpInflationGraph"].as<string>());
}
