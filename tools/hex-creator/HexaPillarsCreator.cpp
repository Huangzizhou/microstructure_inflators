#include "hexlib.h"

void usage(int exitVal) {
    cout << "Usage: HexaPillarsCreator p1 p2 p3 p4 out.wire out.msh" << endl;
    exit(exitVal);
}

int main(int argc, char *argv[])
{
    // Declare the supported options.
    po::options_description desc("Hexa pillars creator");
    desc.add_options()
            ("help", "produce help message")
            ("triangle-side-factor", po::value<double>(), "triangle side factor")
            ("number-pillars", po::value<unsigned>(), "number of pillars")
            ("pillar-area-factor", po::value<double>(), "pillar area factor")
            ("thickness-factor", po::value<double>(), "thickness factor")
            ("output-wire", po::value<std::string>(), "output wire")
            ("output-mesh", po::value<std::string>(), "output msh")
            ;

    po::positional_options_description p;
    p.add("triangle-side-factor", 1);
    p.add("number-pillars", 1);
    p.add("pillar-area-factor", 1);
    p.add("thickness-factor", 1);
    p.add("output-wire", 1);
    p.add("output-mesh", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || vm.size() != 6) {
        usage(1);
        return 1;
    }

    //          . b
    //     . '  |
    // a._______| c
    //       d
    double triangle_side_factor = vm["triangle-side-factor"].as<double>();
    unsigned num_pillars = vm["number-pillars"].as<unsigned>();
    double pillar_area_ratio = vm["pillar-area-factor"].as<double>();
    double thickness_ratio = vm["thickness-factor"].as<double>();
    string out_wire = vm["output-wire"].as<string>();
    string out_mesh = vm["output-mesh"].as<string>();

    double s = 2 / sqrt(3);
    double triangle_side = triangle_side_factor * s * sqrt(3);
    double thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars);

    Matrix<double, 2, Dynamic> vertices;
    vector<vector<int>> edges;
    vector<pair<vector<Point>, double> > custom_pairs;

    cout << "Constructing " + out_wire + " ..." << endl;
    generate_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio,
            vertices, edges, custom_pairs);

    // print wire output
    create_wire(vertices, edges, out_wire);

    cout << "Inflating ..." << endl;

    // Computing void thickness and necessary resolution
    double thickness_void = (triangle_side*pillar_area_ratio - num_pillars*thickness) / (num_pillars - 1);
    int min_resolution = max(2 / thickness_void, 2 / thickness);
    int chosen_resolution = pow(2, ceil(log(min_resolution) / log(2)));

    if (chosen_resolution > 1024) {
        cout << "Resolution of " << chosen_resolution << "is too big" << endl;
        cout << "Skipping experiment!" << endl;
        return 1;
    }

    if (chosen_resolution < 64)
        chosen_resolution = 64;

    cout << "Thickness void: " << thickness_void << endl;
    cout << "Minimum resolution: " << min_resolution <<  endl;
    cout << "Chosen resolution: " << chosen_resolution << endl;

    inflate_hexagonal_box(out_wire, 0.00001, 0.00001, out_mesh, custom_pairs, chosen_resolution);

    return 0;
}

