#include "HexaPillarsInflator.hh"

using namespace std;

HexaPillarsInflator::HexaPillarsInflator(const std::vector<Real> &initial_params, double p2, char structure_type) {
    using Point = Vector2d;

    HexLib<double> hexlib;
    Matrix<double, 2, Dynamic> vertices;
    vector<vector<int>> edges;
    vector<pair<vector<Point>, double> > custom_pairs;
    string out_wire = "temp-hexa-pillars-inflator.wire";

    double triangle_side_factor, pillar_area_ratio, thickness_ratio;
    double center_thickness_ratio, ninja_factor; // new parameters
    double joint_thickness_factor, joint_offset_factor; // newest parameters
    unsigned num_pillars = p2;

    m_structure_type = structure_type;
    if (m_structure_type == '+') {
        triangle_side_factor = initial_params[0];
        pillar_area_ratio = 1.0; // always 1.0 for positive poisson ratios
        thickness_ratio = initial_params[1];

        cout << "Constructing " + out_wire + " ..." << endl;
        hexlib.generate_simpler_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio,
                                                            vertices, edges, custom_pairs);

        m_p1 = triangle_side_factor;
        m_p2 = num_pillars;
        m_p3 = pillar_area_ratio;
        m_p4 = thickness_ratio;
    }
    else {
        triangle_side_factor = initial_params[0];
        pillar_area_ratio = initial_params[1];
        thickness_ratio = initial_params[2];
        center_thickness_ratio = initial_params[3];
        ninja_factor = initial_params[4];
        joint_thickness_factor = initial_params[5];
        joint_offset_factor = initial_params[6];

        cout << "p1: " << triangle_side_factor << endl;
        cout << "p2: " << num_pillars << endl;
        cout << "p3: " << pillar_area_ratio <<  endl;
        cout << "p4: " << thickness_ratio << endl;
        cout << "p5: " << center_thickness_ratio <<  endl;
        cout << "p6: " << ninja_factor << endl;
        cout << "p7: " << joint_thickness_factor <<  endl;
        cout << "p8: " << joint_offset_factor << endl;

        cout << "Constructing " + out_wire + " ..." << endl;
        hexlib.generate_auxetic_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio,
                                                            center_thickness_ratio, ninja_factor, joint_thickness_factor,
                                                            joint_offset_factor, vertices, edges, custom_pairs);

        m_p1 = triangle_side_factor;
        m_p2 = num_pillars;
        m_p3 = pillar_area_ratio;
        m_p4 = thickness_ratio;
        m_p5 = center_thickness_ratio;
        m_p6 = ninja_factor;
        m_p7 = joint_thickness_factor;
        m_p8 = joint_offset_factor;
    }


    hexlib.create_wire(vertices, edges, out_wire);

    // Create and save inflator
    m_infl = Future::make_unique<IsoinflatorWrapper<2>>(out_wire, "doubly_periodic", true, 2);
}

void HexaPillarsInflator::configureResolution(const std::vector<Real> &params) {
    HexLib<double> hexlib;
    double triangle_side_factor;
    unsigned num_pillars;
    double pillar_area_ratio;
    double thickness_ratio;
    double center_thickness_ratio;
    double ninja_factor;
    double joint_thickness_factor, joint_offset_factor; // newest parameters

    double triangle_side, thickness;
    double thickness_void;

    int chosen_resolution;

    if (m_structure_type == '+') {
        triangle_side_factor = params[0];
        num_pillars = m_p2;
        pillar_area_ratio = m_p3;
        thickness_ratio = params[1];

        double s = 2 / sqrt(3);
        triangle_side = triangle_side_factor * s * sqrt(3);
        thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars);
        thickness_void = (triangle_side*pillar_area_ratio - num_pillars*thickness) / (num_pillars - 1);


        cout << "p1: " << triangle_side_factor << endl;
        cout << "p2: " << num_pillars << endl;
        cout << "p3: " << pillar_area_ratio <<  endl;
        cout << "p4: " << thickness_ratio << endl;

        // Computing void thickness and necessary resolution
        double multiplier =  1.0;
        int min_resolution = multiplier * max(2 / thickness_void, 2 / thickness);
        chosen_resolution = pow(2, ceil(log(min_resolution) / log(2)));

        cout << "Multiplier: " << multiplier << endl;
        cout << "Thickness void: " << thickness_void << endl;
        cout << "Minimum resolution: " << min_resolution <<  endl;
    }
    else {
        triangle_side_factor = params[0];
        num_pillars = m_p2;
        pillar_area_ratio = params[1];
        thickness_ratio = params[2];
        center_thickness_ratio = params[3];
        ninja_factor = params[4];
        joint_thickness_factor = params[5];
        joint_offset_factor = params[6];

        double parallelogram_side = 3.0;
        double s = parallelogram_side / 3.0;
        double pillar_area = hexlib.get_pillar_area(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio, center_thickness_ratio, ninja_factor);
        triangle_side = triangle_side_factor * s;
        thickness = hexlib.get_thickness(min(thickness_ratio, center_thickness_ratio), max(thickness_ratio, center_thickness_ratio), 1, num_pillars, pillar_area);
        thickness_void = hexlib.get_spacing(thickness_ratio, center_thickness_ratio, num_pillars, pillar_area);
        double thickness_ninja_void = triangle_side * pillar_area_ratio * ninja_factor; // trying to give some importance to the size of the whole between last pillar and beginning of neighbor triangle
        double thickness_joint_void = hexlib.get_joint_void(triangle_side_factor, num_pillars, pillar_area_ratio, thickness_ratio, center_thickness_ratio, ninja_factor, joint_thickness_factor, joint_offset_factor);

        cout << "p1: " << triangle_side_factor << endl;
        cout << "p2: " << num_pillars << endl;
        cout << "p3: " << pillar_area_ratio <<  endl;
        cout << "p4: " << thickness_ratio << endl;
        cout << "p5: " << center_thickness_ratio <<  endl;
        cout << "p6: " << ninja_factor << endl;
        cout << "p7: " << joint_thickness_factor <<  endl;
        cout << "p8: " << joint_offset_factor << endl;

        // Computing void thickness and necessary resolution
        double multiplier =  3/max(0.5, pillar_area_ratio);
        int min_resolution = multiplier * max(max(max(2 / thickness_void, 2 / thickness), 2 / thickness_ninja_void), 2 / thickness_joint_void);
        chosen_resolution = pow(2, ceil(log(min_resolution) / log(2)));

        cout << "Multiplier: " << multiplier << endl;
        cout << "Thickness void: " << thickness_void << endl;
        cout << "Minimum resolution: " << min_resolution <<  endl;
    }

    if (chosen_resolution > 2048)
        chosen_resolution = 2048;

    if (chosen_resolution < 256)
        chosen_resolution = 256;

    cout << "Chosen resolution: " << chosen_resolution << endl;

    unsigned coarsening = 2;
    unsigned grid_size = pow(2,coarsening) * chosen_resolution;

    // Set values to inflator
    m_infl->meshingOptions().marchingSquaresCoarsening = coarsening;
    m_infl->meshingOptions().marchingSquaresGridSize = grid_size;
    m_infl->meshingOptions().forceMSGridSize = true;
}

void HexaPillarsInflator::m_inflate(const std::vector<Real> &params) {
    cout << "Entering m_inflate" << endl;
    vector<MeshIO::IOVertex> vertices = m_infl->vertices();

    vector<Point2D> vertices_vector;
    for (unsigned i=0; i<vertices.size(); i++) {
        Point2D point2D = vertices[i].point.head(2);
        vertices_vector.push_back(point2D);
    }

    configureResolution(params);

    // Inflate
    m_infl->inflate(hexaPillarsToFullParameters(params));

    if (m_structure_type == '+') {
        m_p1 = params[0];
        m_p4 = params[1];
    }
    else {
        m_p1 = params[0];
        m_p3 = params[1];
        m_p4 = params[2];
        m_p5 = params[3];
        m_p6 = params[4];
        m_p7 = params[5];
        m_p8 = params[6];
    }

    cout << "Exiting m_inflate" << endl;
}


bool HexaPillarsInflator::isPrintable(const std::vector<Real> &params) {
    return m_infl->isPrintable(hexaPillarsToFullParameters(params));
}


std::vector<VectorField<Real, 2>> HexaPillarsInflator::volumeShapeVelocities() const {
    std::vector<VectorField<Real, 2>> original_velocities = m_infl->volumeShapeVelocities();
    const size_t number_hexa_params = numParameters();
    std::vector<VectorField<Real, 2>> result(number_hexa_params);

    if (m_structure_type == '+') {
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

        for (unsigned index_original = 0; index_original < resulting_parameters.size(); index_original++) {
            ADScalar original_param = resulting_parameters[index_original];
            Matrix<double, 2, 1, 0, 2, 1> original_param_derivatives = original_param.derivatives();

            for (unsigned index_hexa = 0; index_hexa < number_hexa_params; index_hexa++) {
                partials[index_hexa][index_original] = original_param_derivatives(index_hexa);
            }
        }

        cout << "Computing velocities by chain rule: " << endl;

        for (size_t index_hexa = 0; index_hexa < number_hexa_params; index_hexa++) {
            auto &vvel = result[index_hexa]; // each element corresponds to the derivatives in each vertex wrt a fixed p_i (i = 1 or 4)
            vvel.resizeDomain(vertices().size());
            vvel.clear();

            vector<double> partials_hexa_param = partials[index_hexa];

            for (unsigned index_original = 0; index_original < resulting_parameters.size(); index_original++) {
                vvel += original_velocities[index_original] * partials_hexa_param[index_original];
            }
        }
    }
    else {
        using PVec = Eigen::Matrix<Real, 7, 1>;
        using ADScalar = Eigen::AutoDiffScalar<PVec>;

        ADScalar p1(m_p1, 7, 0);
        ADScalar p3(m_p3, 7, 1);
        ADScalar p4(m_p4, 7, 2);
        ADScalar p5(m_p5, 7, 3);
        ADScalar p6(m_p6, 7, 4);
        ADScalar p7(m_p7, 7, 5);
        ADScalar p8(m_p8, 7, 6);

        vector<ADScalar> input_params = {p1, p3, p4, p5, p6, p7, p8};
        vector<ADScalar> resulting_parameters =  hexaPillarsToFullParameters<ADScalar>(input_params);

        cout << "Creating partials: " << endl;

        vector<vector<double>> partials(7);

        partials[0].resize(resulting_parameters.size());
        partials[1].resize(resulting_parameters.size());
        partials[2].resize(resulting_parameters.size());
        partials[3].resize(resulting_parameters.size());
        partials[4].resize(resulting_parameters.size());
        partials[5].resize(resulting_parameters.size());
        partials[6].resize(resulting_parameters.size());

        for (unsigned index_original = 0; index_original < resulting_parameters.size(); index_original++) {
            ADScalar original_param = resulting_parameters[index_original];
            Matrix<double, 7, 1, 0, 7, 1> original_param_derivatives = original_param.derivatives();

            for (unsigned index_hexa = 0; index_hexa < number_hexa_params; index_hexa++) {
                partials[index_hexa][index_original] = original_param_derivatives(index_hexa);
            }
        }

        cout << "Computing velocities by chain rule: " << endl;

        for (size_t index_hexa = 0; index_hexa < number_hexa_params; index_hexa++) {
            auto &vvel = result[index_hexa]; // each element corresponds to the derivatives in each vertex wrt a fixed p_i
            vvel.resizeDomain(vertices().size());
            vvel.clear();

            vector<double> partials_hexa_param = partials[index_hexa];

            for (unsigned index_original = 0; index_original < resulting_parameters.size(); index_original++) {
                vvel += original_velocities[index_original] * partials_hexa_param[index_original];
            }
        }

        cout << "Exiting volume shape velocities " << endl;
    }

    return result;
}

template<typename TReal>
std::vector<TReal> HexaPillarsInflator::hexaPillarsToFullParameters(const std::vector<TReal> &hexaPillarsParameters) const {
    cout << "Entering hexaPillarsToFullParameters" << endl;

    Matrix<TReal, 2, Dynamic> vertices;
    vector<vector<int>> edges;
    vector<pair<vector<Eigen::Matrix<TReal, 2, 1>>, TReal> > custom_pairs;
    HexLib<TReal> hexlib;


    if (m_structure_type == '+') {
        TReal triangle_side_factor = hexaPillarsParameters[0];
        unsigned num_pillars = m_p2;
        double pillar_area_ratio = 1.0;
        TReal thickness_ratio = hexaPillarsParameters[1];

        hexlib.generate_simpler_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio,
                                                            thickness_ratio,
                                                            vertices, edges, custom_pairs);
    }
    else {
        TReal triangle_side_factor = hexaPillarsParameters[0];
        unsigned num_pillars = m_p2;
        TReal pillar_area_ratio = hexaPillarsParameters[1];
        TReal thickness_ratio = hexaPillarsParameters[2];
        TReal center_thickness_ratio = hexaPillarsParameters[3];
        TReal ninja_factor = hexaPillarsParameters[4];
        TReal joint_thickness_factor = hexaPillarsParameters[5];
        TReal joint_offset_factor = hexaPillarsParameters[6];

        hexlib.generate_auxetic_topology_and_thickness_info(triangle_side_factor, num_pillars, pillar_area_ratio,
                                                            thickness_ratio, center_thickness_ratio, ninja_factor,
                                                            joint_thickness_factor, joint_offset_factor,
                                                            vertices, edges, custom_pairs);
    }

    vector<Eigen::Matrix<TReal, 2, 1>> vertices_vector;
    for (unsigned i=0; i<vertices.cols(); i++) {
        vertices_vector.push_back(vertices.col(i));
    }

    vector<Eigen::Matrix<TReal, 2, 1>> independent_vertices = hexlib.extract_independent_vertices(vertices_vector);
    vector<TReal> parameters = hexlib.generate_parameters(independent_vertices, 1e-5, 0.0, custom_pairs);

    cout << "Exiting hexaPillarsToFullParameters" << endl;

    return parameters;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
HexaPillarsInflator::selfSupportingConstraints(const std::vector<double> &/*params*/) const {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> emptyMatrix;
    return emptyMatrix;
}

