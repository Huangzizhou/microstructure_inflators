#include <isosurface_inflator/WireMesh.hh>
#include <isosurface_inflator/PatternSignedDistance.hh>
#include <isosurface_inflator/IsosurfaceInflatorConfig.hh>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <openvdb/openvdb.h>
#include <openvdb/tools/SignedFloodFill.h>

#include <iostream>
#include <functional>

namespace po = boost::program_options;
using namespace std;

template<class WMesh>
void execute(const string &wirePath, string paramString, const string &igraphPath, const string &outPath, const int resolution) {
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

    const double bg_val = 5.0 / resolution;
    openvdb::initialize();
    openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(bg_val);
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
    
    const auto bbox = psd.boundingBox();
    Point3d minCorner = bbox.minCorner;
    Point3d maxCorner = bbox.maxCorner;
    std::cout << "minCorner " << minCorner.transpose() << "\n";
    std::cout << "maxCorner " << maxCorner.transpose() << "\n";
    Eigen::Vector3i res(resolution, resolution, resolution);
    Point3d h = (maxCorner - minCorner).array() / res.cast<double>().array();
    int compressed = 0;
    Point3d evalPt;
    openvdb::Coord ijk;
    int &i = ijk[0], &j = ijk[1], &k = ijk[2];
    for (i = 0; i <= res(0); i++)
    {
        evalPt(0) = minCorner(0) + i * h(0);
        for (j = 0; j <= res(1); j++)
        {
            evalPt(1) = minCorner(1) + j * h(1);
            for (k = 0; k <= res(2); k++)
            {
                evalPt(2) = minCorner(2) + k * h(2);

                double val = psd.signedDistance(evalPt);
                if (abs(val) < bg_val)
                    accessor.setValue(ijk, val);
                else
                    compressed++;
            }
        }
    }
    openvdb::tools::signedFloodFill(grid->tree());

    std::cout << "Saved space " << compressed / (double) res.prod() * 100 << " %\n";

    grid->setName("level_set");
    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    openvdb::io::File file(outPath);

    openvdb::GridPtrVec(grids);
    grids.push_back(grid);
    file.write(grids);
    file.close();
}

map<string, function<void(const string &, string, const string &, const string &, const int)>> execImpls =
    {{          "cubic", execute<WireMesh<Symmetry::Cubic<>>         >},
     {    "orthotropic", execute<WireMesh<Symmetry::Orthotropic<>>   >},
     {"triply_periodic", execute<WireMesh<Symmetry::TriplyPeriodic<>>>}};

void usage(int status, const po::options_description &visible_opts) {
    cerr << "Usage: ./SDF_cli pattern_symmetry graph outfile -p parameters [options]" << endl;
    cerr << "eg: ./SDF_cli cubic pattern0746.wire out.vdb -p \"0.25 0.5 0.25 0.25 0.5 0.25 0.25 0.25 0.25\"" << endl;
    cout << visible_opts << endl;
    exit(status);
}

po::variables_map parseCmdLine(int argc, char *argv[]) {
    po::options_description hidden_opts("Hidden Arguments");
    hidden_opts.add_options()
        ("symmetry",  po::value<string>(), "Pattern symmetry (triply_periodic,  orthotropic,  cubic)")
        ("wire"    ,  po::value<string>(), "input wire file")
        ("vdb"    ,  po::value<string>(), "output vdb file")
        ;
    po::positional_options_description p;
    p.add("symmetry", 1);
    p.add("wire", 1);
    p.add("vdb", 1);

    // Options visible in the help message.
    po::options_description visible_opts;
    visible_opts.add_options()("help,h", "Produce this help message")
                              ("dumpInflationGraph,D", po::value<string>()->default_value(""),  " Output the inflation graph to the specified path")
                              ("params,p", po::value<string>(),              " Pattern parameters")
                              ("resolution,r", po::value<int>(),              " SDF resolution")
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

    int resolution = 100;
    if (args.count("resolution")) resolution = args["resolution"].as<int>();

    exec(args["wire"].as<string>(), paramString,
         args["dumpInflationGraph"].as<string>(), args["vdb"].as<string>(), resolution);
}
