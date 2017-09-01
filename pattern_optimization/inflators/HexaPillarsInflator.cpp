#include "HexaPillarsInflator.h"
#include "../../tools/hex-creator/hexlib.h"

using namespace std;

HexaPillarsInflator::HexaPillarsInflator(const std::vector<Real> &initial_params, double p2) {
    cout << "HexaPillarsInflator - enter" << endl;

    Matrix<double, 2, Dynamic> vertices;
    vector<vector<int>> edges;
    vector<pair<vector<Point>, double> > custom_pairs;
    string out_wire = "temp-hexa-pillars-inflator.wire";

    double triangle_side_factor = initial_params[0];
    unsigned num_pillars = p2;
    double pillar_area_ratio = 1.0;
    double thickness_ratio = initial_params[1];

    cout << "Constructing " + out_wire + " ..." << endl;
    generate_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio,
                                         vertices, edges, custom_pairs);


    create_wire(vertices, edges, out_wire);

    // Create and save inflator
    m_infl = Future::make_unique<IsoinflatorWrapper<2>>(out_wire, "doubly_periodic", true, 2);


    m_p1 = triangle_side_factor
    cout << "HexaPillarsInflator - exit" << endl;
}

void HexaPillarsInflator::configureResolution(const std::vector<Real> &params) {
    double triangle_side_factor = params[0];
    unsigned num_pillars = m_p2;
    double pillar_area_ratio = m_p3;
    double thickness_ratio = params[1];

    double parallelogram_side = 2.0;
    double l = parallelogram_side / 2.0;
    double s = 2 / sqrt(3);
    double triangle_side = triangle_side_factor * s * sqrt(3);
    double thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars);

    // Computing void thickness and necessary resolution
    double thickness_void = (triangle_side*pillar_area_ratio - num_pillars*thickness) / (num_pillars - 1);
    int min_resolution = max(2 / thickness_void, 2 / thickness);
    int chosen_resolution = pow(2, ceil(log(min_resolution) / log(2)));

    if (chosen_resolution > 1024)
        chosen_resolution = 1024;

    if (chosen_resolution < 64)
        chosen_resolution = 64;

    cout << "Thickness void: " << thickness_void << endl;
    cout << "Minimum resolution: " << min_resolution <<  endl;
    cout << "Chosen resolution: " << chosen_resolution << endl;

    unsigned coarsening = 2;
    unsigned grid_size = pow(2,coarsening) * chosen_resolution;

    // Set values to inflator
    m_infl->meshingOptions().marchingSquaresCoarsening = coarsening;
    m_infl->meshingOptions().marchingSquaresGridSize = grid_size;
    m_infl->meshingOptions().forceMSGridSize = true;

}

void HexaPillarsInflator::m_inflate(const std::vector<Real> &params) {
    vector<MeshIO::IOVertex> vertices = m_infl->vertices();

    vector<Point2D> vertices_vector;
    for (unsigned i=0; i<vertices.size(); i++) {
        Point3D point3D = vertices[i].point;
        Point2D point2D = vertices[i].point.head(2);
        vertices_vector.push_back(point2D);
    }

    configureResolution(params);

    // Inflate
    m_infl->inflate(hexaPillarsToFullParameters(params));

    m_p1 = params[0];
    m_p3 = params[1];
}


bool HexaPillarsInflator::isPrintable(const std::vector<Real> &params) {
    return m_infl->isPrintable(hexaPillarsToFullParameters(params));
}


std::vector<VectorField<Real, 2>> HexaPillarsInflator::volumeShapeVelocities() const {
    //TODO: compute the derivatives of the old parameters with respect to the new ones.
    // And use chain rule!

    using PVec = Eigen::Matrix<Real, 3, 1>;
    using ADScalar = Eigen::AutoDiffScalar<PVec>;

    ADScalar p1(m_p1, 2, 0);
    ADScalar p3(m_p3, 2, 1);


    for (int index=0; index<m_params.size(); index++) {
        Real cur_param = m_params[index];



    }


    return m_infl->volumeShapeVelocities();
}

template<typename T>
std::vector<T> HexaPillarsInflator::hexaPillarsToFullParameters(const std::vector<T> &hexaPillarsParameters) {
    Matrix<double, 2, Dynamic> vertices;
    vector<vector<int>> edges;
    vector<pair<vector<Point>, double> > custom_pairs;

    cout << "Entering hexaPillarsToFullParameters" << endl;

    double triangle_side_factor = hexaPillarsParameters[0];
    unsigned num_pillars = m_params[1]; //TODO: static_cast<unsigned>(hexaPillarsParameters[1]);
    double pillar_area_ratio = 1.0; //TODO: hexaPillarsParameters[2];
    double thickness_ratio = hexaPillarsParameters[3];

    generate_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio,
                                         vertices, edges, custom_pairs);

    vector<Point> vertices_vector;
    for (unsigned i=0; i<vertices.cols(); i++) {
        vertices_vector.push_back(vertices.col(i));
    }

    //Check if vertices are the same as in the inflator
    vector<MeshIO::IOVertex> inflator_vertices = m_infl->vertices();
    for (unsigned i=0; i<inflator_vertices.size(); i++) {
        Point3D inflator_point = inflator_vertices[i].point;
        Point2D new_point = vertices_vector[i];

        assert(abs(inflator_point(0) - new_point(0)) < 1e-7);
        assert(abs(inflator_point(1) - new_point(1)) < 1e-7);
    }

    vector<Point> independent_vertices = extract_independent_vertices(vertices_vector);
    vector<double> parameters = generate_parameters(independent_vertices, 1e-5, 1e-5, custom_pairs);

    cout << "Exiting hexaPillarsToFullParameters" << endl;

    return parameters;
}

