#include "HexaPillarsInflator.h"

using namespace std;

HexaPillarsInflator::HexaPillarsInflator(const std::vector<Real> &initial_params, double p2) {
    using Point = Vector2d;

    HexLib<double> hexlib;
    Matrix<double, 2, Dynamic> vertices;
    vector<vector<int>> edges;
    vector<pair<vector<Point>, double> > custom_pairs;
    string out_wire = "temp-hexa-pillars-inflator.wire";

    double triangle_side_factor = initial_params[0];
    unsigned num_pillars = p2;
    double pillar_area_ratio = 1.0; // always 1.0 for positive poisson ratios
    double thickness_ratio = initial_params[1];

    cout << "Constructing " + out_wire + " ..." << endl;
    hexlib.generate_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio,
                                        vertices, edges, custom_pairs);

    hexlib.create_wire(vertices, edges, out_wire);

    // Create and save inflator
    m_infl = Future::make_unique<IsoinflatorWrapper<2>>(out_wire, "doubly_periodic", true, 2);

    m_p1 = triangle_side_factor;
    m_p2 = num_pillars;
    m_p3 = pillar_area_ratio;
    m_p4 = thickness_ratio;
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

    if (chosen_resolution < 32)
        chosen_resolution = 32;

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
    m_p4 = params[1];
}


bool HexaPillarsInflator::isPrintable(const std::vector<Real> &params) {
    return m_infl->isPrintable(hexaPillarsToFullParameters(params));
}


std::vector<VectorField<Real, 2>> HexaPillarsInflator::volumeShapeVelocities() const {
    std::vector<VectorField<Real, 2>> original_velocities = m_infl->volumeShapeVelocities();
    const size_t number_hexa_params = numParameters();
    std::vector<VectorField<Real, 2>> result(number_hexa_params);

    using PVec = Eigen::Matrix<Real, 2, 1>;
    using ADScalar = Eigen::AutoDiffScalar<PVec>;

    ADScalar p1(m_p1, 2, 0);
    ADScalar p4(m_p4, 2, 1);

    vector<ADScalar> input_params = {p1, p4};
    vector<ADScalar> resulting_parameters =  hexaPillarsToFullParameters<ADScalar>(input_params);

    cout << "Creating partials: " << endl;

    vector<vector<double>> partials(2);
    partials[0].resize(resulting_parameters.size());
    partials[1].resize(resulting_parameters.size());

    for (int index_original = 0; index_original < resulting_parameters.size(); index_original++) {
        ADScalar original_param = resulting_parameters[index_original];
        Matrix<double, 2, 1, 0, 2, 1> original_param_derivatives = original_param.derivatives();

        for (int index_hexa = 0; index_hexa < number_hexa_params; index_hexa++) {
            partials[index_hexa][index_original] = original_param_derivatives(index_hexa);
        }
    }

    cout << "Computing velocities by chain rule: " << endl;

    for (size_t index_hexa = 0; index_hexa < number_hexa_params; index_hexa++) {
        auto &vvel = result[index_hexa]; // each element corresponds to the derivatives in each vertex wrt a fixed p_i (i = 1 or 4)
        vvel.resizeDomain(vertices().size());
        vvel.clear();

        vector<double> partials_hexa_param = partials[index_hexa];

        for (int index_original = 0; index_original < resulting_parameters.size(); index_original++) {
            vvel += original_velocities[index_original] * partials_hexa_param[index_original];
        }
    }

    return result;
}

template<typename TReal>
std::vector<TReal> HexaPillarsInflator::hexaPillarsToFullParameters(const std::vector<TReal> &hexaPillarsParameters) const {
    Matrix<TReal, 2, Dynamic> vertices;
    vector<vector<int>> edges;
    vector<pair<vector<Eigen::Matrix<TReal, 2, 1>>, TReal> > custom_pairs;
    HexLib<TReal> hexlib;

    TReal triangle_side_factor = hexaPillarsParameters[0];
    unsigned num_pillars = m_p2; //TODO: static_cast<unsigned>(hexaPillarsParameters[1]);
    double pillar_area_ratio = 1.0; //TODO: hexaPillarsParameters[2];
    TReal thickness_ratio = hexaPillarsParameters[1];

    hexlib.generate_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio,
                                         vertices, edges, custom_pairs);

    vector<Eigen::Matrix<TReal, 2, 1>> vertices_vector;
    for (unsigned i=0; i<vertices.cols(); i++) {
        vertices_vector.push_back(vertices.col(i));
    }

    vector<Eigen::Matrix<TReal, 2, 1>> independent_vertices = hexlib.extract_independent_vertices(vertices_vector);
    vector<TReal> parameters = hexlib.generate_parameters(independent_vertices, 1e-5, 1e-5, custom_pairs);

    return parameters;
}

