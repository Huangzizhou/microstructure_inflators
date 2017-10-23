#include "hexlib.h"

typedef Vector2d APoint;

void usage(int exitVal) {
    cout << "Usage: AuxeticHexaPillarsCreator p1 p2 p3 p4 p5 p6 p7 p8 out.wire out.msh" << endl;
    exit(exitVal);
}


int main(int argc, char *argv[]) {
    HexLib<double> hexlib;
    Matrix<double, 2, Dynamic> vertices;
    vector<vector<int>> edges;
    vector<pair<vector<APoint>, double> > custom_pairs;

    // Declare the supported options.
    po::options_description desc("Auxetic hexa pillars creator");
    desc.add_options()
            ("help", "produce help message")
            ("triangle-side-factor", po::value<double>(), "triangle side factor")
            ("number-pillars", po::value<unsigned>(), "number of pillars")
            ("pillar-area-factor", po::value<double>(), "pillar area factor")
            ("min-thickness-factor", po::value<double>(), "min thickness factor")
            ("max-thickness-factor", po::value<double>(), "max thickness factor")
            ("ninja-factor", po::value<double>(), "ninja factor")
            ("joint-thickness-factor", po::value<double>(), "joint thickness factor")
            ("joint-offset-factor", po::value<double>(), "joint offset factor")
            ("output-wire", po::value<std::string>(), "output wire")
            ("output-mesh", po::value<std::string>(), "output msh");

    po::positional_options_description p;
    p.add("triangle-side-factor", 1);
    p.add("number-pillars", 1);
    p.add("pillar-area-factor", 1);
    p.add("min-thickness-factor", 1);
    p.add("max-thickness-factor", 1);
    p.add("ninja-factor", 1);
    p.add("joint-thickness-factor", 1);
    p.add("joint-offset-factor", 1);
    p.add("output-wire", 1);
    p.add("output-mesh", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") || vm.size() != 10) {
        usage(1);
        return 1;
    }


    //       .b
    //      /\
    //     /  \
    //    /    \
    // a /______\ c
    //
    double triangle_side_ratio = vm["triangle-side-factor"].as<double>();
    unsigned num_pillars = vm["number-pillars"].as<unsigned>();
    double pillar_area_ratio = vm["pillar-area-factor"].as<double>();
    double min_thickness_ratio = vm["min-thickness-factor"].as<double>();
    double max_thickness_ratio = vm["max-thickness-factor"].as<double>();
    double ninja_factor = vm["ninja-factor"].as<double>();
    double joint_thickness_factor = vm["joint-thickness-factor"].as<double>();
    double joint_offset_factor = vm["joint-offset-factor"].as<double>();
    string out_wire = vm["output-wire"].as<string>();
    string out_mesh = vm["output-mesh"].as<string>();

    double p1 = triangle_side_ratio;
    double p2 = num_pillars;
    double p3 = pillar_area_ratio;
    double p4 = min_thickness_ratio;
    double p5 = max_thickness_ratio;
    double p6 = ninja_factor;
    double p7 = joint_thickness_factor;
    double p8 = joint_offset_factor;

    double parallelogram_side = 3.0;
    double s = parallelogram_side / 3.0;
    double triangle_side = triangle_side_ratio * s;
    double triangle_y_position = s * sqrt(3)/2.0 * (1-triangle_side_ratio);

    // define important vertices of simplex used to build entire parallelogram structure
    Eigen::Matrix<double, 2, 1> origin, a, b, c;
    origin << 0, 0;
    a << 0, 0;
    b << s/2.0, s * sqrt(3)/2.0;
    c << s, 0;

    // define vertices of triangle
    Eigen::Matrix<double, 2, 1> q1, q2, w, z;
    Eigen::Matrix<double, 2, 1> q1_reflected, q2_reflected, w_reflected, z_reflected;
    Eigen::Matrix<double, 2, 1> ba_unit, ab_unit;
    q1 << s/2.0 * (1-triangle_side_ratio), triangle_y_position;
    q2 << s/2.0 * (1+triangle_side_ratio), triangle_y_position;

    w  << q2[0] - p1*p3*s, triangle_y_position;
    ba_unit = (b - a) / (b-a).norm();
    z = ba_unit * triangle_side_ratio*p3*s*(1-p6) + w;

    q1_reflected << q2[0], -q2[1];
    q2_reflected << q1[0], -q1[1];
    w_reflected  << q2_reflected[0] + p1*p3*s , -triangle_y_position;
    ab_unit = (a-b) / (a-b).norm();
    z_reflected = ab_unit * triangle_side_ratio*p3*s*(1-p6) + w_reflected;

    hexlib.generate_auxetic_topology_and_thickness_info(triangle_side_ratio, num_pillars, pillar_area_ratio, min_thickness_ratio, max_thickness_ratio, ninja_factor, joint_thickness_factor, joint_offset_factor, vertices, edges, custom_pairs);

    // FINALLY, CREATE WIRE
    hexlib.create_wire(vertices, edges, out_wire);

    // computing thickness and spacing
    double pillar_area = (z - q2).norm();
    double new_pillar_area = hexlib.get_pillar_area(triangle_side_ratio, num_pillars, pillar_area_ratio, min_thickness_ratio, max_thickness_ratio, ninja_factor);
    if (abs(pillar_area - new_pillar_area) > TOLERANCE) {
        cout << "WRONG COMPUTATION OF PILLAR AREA?" << endl;
        exit(1);
    }
    assert(abs(pillar_area - new_pillar_area) < TOLERANCE);
    double thickness = hexlib.get_thickness(min(min_thickness_ratio, max_thickness_ratio), max(min_thickness_ratio, max_thickness_ratio), 1, num_pillars, pillar_area);
    double spacing = hexlib.get_spacing(min_thickness_ratio, max_thickness_ratio, num_pillars, pillar_area);

    double thickness_void = spacing;
    double min_resolution =  2 / max(0.5, pillar_area_ratio) * max(3 / thickness_void, 3 / thickness);
    double chosen_resolution = pow(2, ceil(log(min_resolution) / log(2)));

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

    hexlib.inflate_hexagonal_box(out_wire, 0.00001, 0.0000, out_mesh, custom_pairs, chosen_resolution);

    return 0;
}
