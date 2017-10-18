//
// Created by Davi Colli Tozoni on 8/29/17.
//

#ifndef MICROSTRUCTURES_HEXLIB_H
#define MICROSTRUCTURES_HEXLIB_H

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <cmath>
#include <boost/program_options.hpp>
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include "json.hpp"

#define TOLERANCE 1e-7

using json = nlohmann::json;
using Eigen::Vector2d;
namespace po = boost::program_options;
using namespace std;
using namespace Eigen;

template <typename TReal>
class HexLib {
public:
    using Point = Eigen::Matrix<TReal, 2, 1>;

    HexLib() {}
    ~HexLib() {}

    std::string execute_command(string cmd) {
        char buffer[128];
        std::string result = "";

        FILE *pipe = popen(cmd.c_str(), "r");
        if (!pipe) throw std::runtime_error("popen() failed!");
        try {
            while (!feof(pipe)) {
                if (fgets(buffer, 128, pipe) != NULL)
                    result += buffer;
            }
        } catch (...) {
            pclose(pipe);
            throw;
        }
        pclose(pipe);

        return result;
    }


    void create_wire(Matrix<double, 2, Dynamic> &vertices, vector<vector<int>> &edges, string out_wire) {
        std::ofstream out_file(out_wire);


        for (unsigned index = 0; index < vertices.cols(); index++) {
            Point vertex = vertices.col(index);
            out_file << "v " << std::fixed << std::setprecision(17) << vertex(0) << " " << vertex(1) << " 0" << endl;
        }

        for (unsigned index = 0; index < edges.size(); index++) {
            vector<int> edge = edges[index];
            out_file << "l " << (edge[0] + 1) << " " << (edge[1] + 1) << endl;
        }

        out_file.close();
    }

    string create_custom_meshing_file(unsigned resolution) {
        string custom_name = "custom_meshing_file.json";
        string original_name = "refined-meshing_opts.json";
        int coarsening = 2;
        json meshing_opts;
        std::ifstream patchFile(original_name);

        // parse
        try {
            patchFile >> meshing_opts;
        } catch (...) {
            std::cerr << "Error parsing the json file" << std::endl;
            return "";
        }

        // change
        unsigned grid_size = pow(2, coarsening) * resolution;
        meshing_opts["marchingSquaresGridSize"] = grid_size;
        meshing_opts["marchingSquaresCoarsening"] = coarsening;

        // write to file
        std::ofstream out_file(custom_name);
        out_file << meshing_opts << std::endl;

        return custom_name;
    }

    Matrix<TReal, 2, Dynamic> reflect(double angle, Matrix<TReal, 2, Dynamic> vertices) {
        double theta = M_PI * (angle / 180);
        Matrix<double, 2, 2> Rl;
        Rl << cos(2 * theta), sin(2 * theta),
                sin(2 * theta), -cos(2 * theta);

        Matrix<TReal, 2, Dynamic> reflected_vertices = Rl * vertices;

        return reflected_vertices;
    }

    Matrix<TReal, 2, Dynamic> rotate(double angle, Matrix<TReal, 2, Dynamic> vertices) {
        double theta = M_PI * (angle / 180);
        Matrix<double, 2, 2> R;
        R << cos(theta), -sin(theta),
                sin(theta), cos(theta);

        Matrix<TReal, 2, Dynamic> rotated_vertices = R * vertices;

        return rotated_vertices;
    }

    Matrix<TReal, 2, Dynamic> push_back(Matrix<TReal, 2, Dynamic> &matrix, Point column) {
        matrix.conservativeResize(Eigen::NoChange, matrix.cols() + 1);

        matrix.col(matrix.cols() - 1) = column;

        return matrix;
    }

    bool equal_vertices(Point v1, Point v2) {
        if ((abs(v1(0) - v2(0)) < TOLERANCE) && (abs(v1(1) - v2(1)) < TOLERANCE))
            return true;

        return false;
    }

    int find_vertex_in_list(vector<Point> &vertices, Point v) {
        int vertex_position = -1;

        for (unsigned int index = 0; index < vertices.size(); index++) {
            Point vertex = vertices[index];
            if (equal_vertices(v, vertex)) {
                vertex_position = index;
                break;
            }
        }
        return vertex_position;
    }

    int find_vertex_in_list(Matrix<TReal, 2, Dynamic> &vertices, Point v) {
        int vertex_position = -1;

        for (unsigned int index = 0; index < vertices.cols(); index++) {
            Point vertex = vertices.col(index);
            if (equal_vertices(v, vertex)) {
                vertex_position = index;
                break;
            }
        }
        return vertex_position;
    }

    int find_edge_in_list(vector<vector<int>> &edges, vector<int> e) {
        int new_v1 = e[0];
        int new_v2 = e[1];

        int edge_position = -1;

        for (unsigned index = 0; index < edges.size(); index++) {
            vector<int> edge = edges[index];

            int v1 = edge[0];
            int v2 = edge[1];

            if ((new_v1 == v1) && (new_v2 == v2)) {
                edge_position = index;
                break;
            }

            if ((new_v1 == v2) && (new_v2 == v1)) {
                edge_position = index;
                break;
            }
        }

        return edge_position;
    }


    TReal min_distance_point_line(Point point, vector<Point> line) {
        TReal x0 = point[0];
        TReal y0 = point[1];
        TReal x1 = line[0](0);
        TReal y1 = line[0](1);
        TReal x2 = line[1](0);
        TReal y2 = line[1](1);

        TReal line_norm = sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));

        TReal distance = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / line_norm;

        return distance;
    }


    Point triangle_incenter(Matrix<TReal, 2, Dynamic> triangle) {
        Point p_a = triangle.col(0);
        Point p_b = triangle.col(1);
        Point p_c = triangle.col(2);

        TReal a = (p_b - p_c).norm();
        TReal b = (p_a - p_c).norm();
        TReal c = (p_a - p_b).norm();

        Point incenter = (a * p_a + b * p_b + c * p_c) / (a + b + c);

        return incenter;
    }

    vector<vector<Point>> polygon_to_edges_descriptions(Matrix<TReal, 2, -1, 0, 2, -1> polygon) {
        vector<vector<Point>> edges_descriptions;

        for (unsigned i = 0; i < polygon.cols(); i++) {
            Point p1 = polygon.col(i);
            Point p2 = polygon.col((i + 1) % polygon.cols());

            edges_descriptions.push_back({p1, p2});
        }

        return edges_descriptions;
    }

    void add_new_edges(vector<vector<Point> > &edges_descriptions,  Matrix<TReal, 2, Dynamic> &current_vertices,
                       vector<vector<int>> &current_edges) {

        for (unsigned index = 0; index < edges_descriptions.size(); index++) {
            vector<Point> edge_description = edges_descriptions[index];
            Point vertex1 = edge_description[0];
            Point vertex2 = edge_description[1];

            int position1 = find_vertex_in_list(current_vertices, vertex1);
            if (position1 == -1) {
                position1 = current_vertices.cols();
                push_back(current_vertices, vertex1);
            }

            int position2 = find_vertex_in_list(current_vertices, vertex2);
            if (position2 == -1) {
                position2 = current_vertices.cols();
                push_back(current_vertices, vertex2);
            }

            vector<int> new_edge = {position1, position2};
            assert(position1 != position2);

            int edge_position = find_edge_in_list(current_edges, new_edge);
            if (edge_position == -1)
                current_edges.push_back(new_edge);
        }
    }

    vector<pair<vector<Point>, TReal> >
    add_polygons_incenters(vector<Matrix<TReal, 2, Dynamic> > polygons,  Matrix<TReal, 2, Dynamic> &vertices,
                           vector<vector<int>> &edges) {
        vector<pair<vector<Point>, TReal> > incenters_thickness_pairs;

        for (unsigned index = 0; index < polygons.size(); index++) {
            Matrix<TReal, 2, Dynamic> pol = polygons[index];

            Point incenter = triangle_incenter(pol);
            TReal thickness = min(min(min_distance_point_line(incenter, {pol.col(0), pol.col(1)}),
                                       min_distance_point_line(incenter, {pol.col(1), pol.col(2)})),
                                   min_distance_point_line(incenter, {pol.col(2), pol.col(0)}));

            incenters_thickness_pairs.push_back({{incenter}, thickness});

            vector<vector<Point>> edges_descriptions;

            for (int i = 0; i < pol.cols(); i++) {
                Point a = pol.col(i);
                Point b = incenter;
                edges_descriptions.push_back({a, b});
            }

            add_new_edges(edges_descriptions, vertices, edges);
        }

        return incenters_thickness_pairs;
    }

    vector<Point> extract_vertices_from_edges_descriptions(vector<vector<Point> > edges_descriptions) {
        vector<Point> extracted_vertices;

        for (unsigned int index = 0; index < edges_descriptions.size(); index++) {
            vector<Point> edge_description = edges_descriptions[index];

            Point vertex1 = edge_description[0];
            Point vertex2 = edge_description[1];

            int position1 = find_vertex_in_list(extracted_vertices, vertex1);
            if (position1 == -1)
                extracted_vertices.push_back(vertex1);

            int position2 = find_vertex_in_list(extracted_vertices, vertex2);
            if (position2 == -1)
                extracted_vertices.push_back(vertex2);
        }

        return extracted_vertices;
    }


    void create_subdivided_segment_with_constant_spacing(Point segment_start, Point segment_end, int n, TReal spacing,
                                                         vector<vector<Point> > &edges_descriptions) {
        Point vector = (segment_end - segment_start);
        Point direction = vector / vector.norm();

        Point last_point = segment_start;

        for (int t = 1; t < (n - 1); t++) {
            Point new_point = segment_start + spacing * t * direction;
            edges_descriptions.push_back({last_point, new_point});
            last_point = new_point;
        }

        edges_descriptions.push_back({last_point, segment_end});
    }

    void add_new_vertices_and_edges(Matrix<TReal, 2, Dynamic> reflected_vertices, vector<vector<int>> reflected_edges,
                                    Matrix<TReal, 2, Dynamic> &vertices, vector<vector<int>> &edges) {
        vector<vector<Point> > new_edges_description;

        for (vector<vector<int>>::iterator it = reflected_edges.begin(); it != reflected_edges.end(); it++) {
            vector<int> edge = *it;

            int v1 = edge[0];
            int v2 = edge[1];

            Point reflected_vertex1 = reflected_vertices.col(v1);
            Point reflected_vertex2 = reflected_vertices.col(v2);

            new_edges_description.push_back({reflected_vertex1, reflected_vertex2});
        }

        add_new_edges(new_edges_description, vertices, edges);
    }

    void simplex_to_whole_parallelogram( Matrix<TReal, 2, Dynamic> &vertices, vector<vector<int>> &edges,
                                        double parallelogram_side) {
        double l = parallelogram_side / 2;

        // Operation 1: reflect against 30 degrees line
        Matrix<TReal, 2, Dynamic> reflected_vertices = reflect(30, vertices);
        vector<vector<int>> reflected_edges;
        std::copy(edges.begin(), edges.end(), std::back_insert_iterator<vector<vector<int> > >(reflected_edges));
        add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges);
        Matrix<TReal, 2, Dynamic> vertices_step_1;
        vector<vector<int>> edges_step_1;
        vertices_step_1 = vertices;
        std::copy(edges.begin(), edges.end(),
                  std::back_insert_iterator<vector<vector<int>>>(edges_step_1));

        // Operation 2: reflect against y axis at x = l
        Point offset_point;
        offset_point << l, 0;
        vertices = vertices.colwise() - offset_point;
        reflected_vertices = reflect(90, vertices);
        reflected_edges.clear();
        std::copy(edges.begin(), edges.end(), std::back_insert_iterator<vector<vector<int> > >(reflected_edges));
        add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges);
        vertices = vertices.colwise() + offset_point;

        // Operation 3: rotate 240 degrees and plug on top of current structure
        offset_point << l, sqrt(3) / 3;
        Matrix<TReal, 2, Dynamic> translated_vertices = vertices_step_1.colwise() - offset_point;
        Matrix<TReal, 2, Dynamic> rotated_vertices = rotate(240, translated_vertices);
        vector<vector<int>> rotated_edges;
        std::copy(edges_step_1.begin(), edges_step_1.end(),
                  std::back_insert_iterator<vector<vector<int> > >(rotated_edges));
        rotated_vertices = rotated_vertices.colwise() + offset_point;
        add_new_vertices_and_edges(rotated_vertices, rotated_edges, vertices, edges);

        // Operation 4:
        offset_point << 2 * l, 0;
        translated_vertices = vertices.colwise() - offset_point;
        reflected_vertices = reflect(120, translated_vertices);
        reflected_edges.clear();
        std::copy(edges.begin(), edges.end(), std::back_insert_iterator<vector<vector<int> > >(reflected_edges));
        reflected_vertices = reflected_vertices.colwise() + offset_point;
        add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges);

        // Operation 5: transform parallelogram to square
        Matrix<TReal, 2, 2> transformation_matrix;
        transformation_matrix << 1.0, -0.577350269189626,
                0.0, 1.154700538379251;
        Matrix<TReal, 2, Dynamic> resulting_vertices = transformation_matrix * vertices;
        vertices = resulting_vertices;

        // Finally, transpose to origin
        offset_point << l, l;
        vertices = vertices.colwise() - offset_point;
    }


    vector<Matrix<TReal, 2, Dynamic> >
    simplex_polygon_to_whole_parallelogram(Matrix<TReal, 2, Dynamic> simplex_polygon, double parallelogram_side) {
        double l = parallelogram_side / 2;
        Matrix<TReal, 2, Dynamic> vertices = simplex_polygon;
        vector<Matrix<TReal, 2, Dynamic> > polygons;
        polygons.emplace_back(simplex_polygon);

        // Operation 1: reflect against 30 degrees line
        Matrix<TReal, 2, Dynamic> reflected_vertices = reflect(30, vertices);
        polygons.emplace_back(reflected_vertices);  // polygons now has 2 pols
        vector<Matrix<TReal, 2, Dynamic> > polygons_step_1;
        std::copy(polygons.begin(), polygons.end(),
                  std::back_insert_iterator<vector<Matrix<TReal, 2, Dynamic> > >(polygons_step_1));

        // Operation 2: reflect against y axis at x = l
        vector<Matrix<TReal, 2, Dynamic> > polygons_copy;
        std::copy(polygons.begin(), polygons.end(),
                  std::back_insert_iterator<vector<Matrix<TReal, 2, Dynamic> > >(polygons_copy));
        for (typename std::vector< Matrix<TReal, 2, Dynamic> >::iterator it = polygons_copy.begin();
             it != polygons_copy.end(); ++it) {
            Point offset_point;
            offset_point << l, 0;
            Matrix<TReal, 2, Dynamic> pol_vertices = (*it).colwise() - offset_point;
            reflected_vertices = reflect(90, pol_vertices);
            reflected_vertices = reflected_vertices.colwise() + offset_point;
            polygons.emplace_back(reflected_vertices); // polygons now has 4 pols
        }

        // Operation 3: rotate 240 degrees and plug on top of current structure
        polygons_copy.clear();
        std::copy(polygons_step_1.begin(), polygons_step_1.end(),
                  std::back_insert_iterator<vector<Matrix<TReal, 2, Dynamic> > >(polygons_copy));
        for (typename std::vector<Matrix<TReal, 2, Dynamic> >::iterator it = polygons_copy.begin();
             it != polygons_copy.end(); ++it) {
            Point offset_point;
            offset_point << l, sqrt(3) / 3;
            Matrix<TReal, 2, Dynamic> translated_vertices = (*it).colwise() - offset_point;
            Matrix<TReal, 2, Dynamic> rotated_vertices = rotate(240, translated_vertices);
            rotated_vertices = rotated_vertices.colwise() + offset_point;
            polygons.emplace_back(rotated_vertices); // polygons now has 6 pols
        }

        // Operation 4:
        polygons_copy.clear();
        std::copy(polygons.begin(), polygons.end(),
                  std::back_insert_iterator<vector<Matrix<TReal, 2, Dynamic> > >(polygons_copy));
        for (typename std::vector<Matrix<TReal, 2, Dynamic> >::iterator it = polygons_copy.begin();
             it != polygons_copy.end(); ++it) {
            Point offset_point;
            offset_point << 2 * l, 0;
            Matrix<TReal, 2, Dynamic> translated_vertices = (*it).colwise() - offset_point;
            reflected_vertices = reflect(120, translated_vertices);
            reflected_vertices = reflected_vertices.colwise() + offset_point;
            polygons.emplace_back(reflected_vertices); // polygons now has 12 pols
        }

        // Operation 5: transform parallelogram to square
        Matrix<TReal, 2, 2> transformation_matrix;
        transformation_matrix << 1.0, -0.577350269189626,
                0.0, 1.154700538379251;
        for (unsigned int index = 0; index < polygons.size(); index++) {
            Matrix<TReal, 2, Dynamic> vertices = polygons[index];
            Matrix<TReal, 2, Dynamic> resulting_vertices = transformation_matrix * vertices;
            Point offset_point;
            offset_point << l, l;
            polygons[index] = resulting_vertices.colwise() - offset_point;
        }

        return polygons;
    }

    bool is_zero(TReal number) {
        if (abs(number) < TOLERANCE)
            return true;

        return false;
    }

    Point independent_vertex_position(Point v) {
        if (is_zero(abs(v(0) - 1.0)))
            v(0) = -1.0;

        if (is_zero(abs(v[1] - 1.0)))
            v(1) = -1.0;

        return v;
    }

    vector<Point> extract_independent_vertices(vector<Point> &vertices) {
        vector<Point> independent_vertices;

        for (unsigned index = 0; index < vertices.size(); index++) {
            Point vertex = vertices[index];
            Point independent_vertex = independent_vertex_position(vertex);

            if (equal_vertices(vertex, independent_vertex))
                independent_vertices.push_back(vertex);
        }

        return independent_vertices;
    }

    vector<TReal> generate_position_parameters(vector<Point> vertices) {
        vector<TReal> parameters;

        for (unsigned i = 0; i < vertices.size(); i++) {
            Point vertex = vertices[i];

            // first, normalize vertex
            Point normalized_vertex;
            normalized_vertex(0) = (vertex(0) + 1.0) / 2;
            normalized_vertex(1) = (vertex(1) + 1.0) / 2;

            // then, verify if each coordinate is free and, if so, compute parameter
            if (!is_zero(normalized_vertex(0)) && !is_zero(normalized_vertex(0) - 1.0)) {
                parameters.push_back(normalized_vertex(0));
            }

            if (!is_zero(normalized_vertex(1)) && !is_zero(normalized_vertex(1) - 1.0)) {
                parameters.push_back(normalized_vertex(1));
            }

            // cause z is always at 0 coordinate
            parameters.push_back(0.5);
        }

        return parameters;
    }

    vector<TReal> generate_parameters(vector<Point> independent_vertices, double default_thickness, double default_blending,
                        vector<pair<vector<Point>, TReal>> custom_thickness_pairs) {
        vector<TReal> parameters = generate_position_parameters(independent_vertices);

        for (unsigned i = 0; i < independent_vertices.size(); i++) {
            parameters.push_back(default_thickness);
        }

        for (unsigned i = 0; i < independent_vertices.size(); i++) {
            parameters.push_back(default_blending);
        }

        size_t offset_thickeness = parameters.size() - 2 * independent_vertices.size();

        // Now, for the extra parameters, find vertices and apply customized thickness
        for (unsigned index = 0; index < custom_thickness_pairs.size(); index++) {
            pair<vector<Point>, TReal> thickness_pair = custom_thickness_pairs[index];

            vector<Point> custom_nodes = thickness_pair.first;
            TReal custom_thickness = thickness_pair.second;

            for (typename vector<Point>::iterator node_it = custom_nodes.begin(); node_it != custom_nodes.end(); node_it++) {
                Point independent_vertex = independent_vertex_position(*node_it);
                int node_index = find_vertex_in_list(independent_vertices, independent_vertex);
                if (node_index >= 0)
                    parameters[offset_thickeness + node_index] = custom_thickness;
                else
                    cout << "[Warning] vertex with customized thickness not found" << endl;
            }
        }

        return parameters;
    }

    void inflate_hexagonal_box(string input_path, double default_thickness, double default_blending, string out_path,
                               vector<pair<vector<Point>, double>> custom_thickness_pairs = vector<pair<vector<Point>, double>>(),
                               unsigned resolution = 64) {

        // Discover the vertices
        vector<Point> vertices;

        std::ifstream input_wire(input_path);

        std::string line;
        while (std::getline(input_wire, line))
        {
            std::istringstream iss(line);
            char prefix;
            double x, y, z;

            // if line does not follow the patttern we expect, just skip it!
            if (!(iss >> prefix >> x >> y >> z)) {
                continue;
            }

            if (prefix == 'v' && abs(z) < TOLERANCE) {
                vertices.push_back({x, y});
            }
        }

        vector<Point> independent_vertices = extract_independent_vertices(vertices);

        vector<double> parameters = generate_parameters(independent_vertices, default_thickness,
                                                                default_blending, custom_thickness_pairs);

        // Write parameters to file
        string parameters_string;
        for (float s : parameters) {
            std::ostringstream buff;
            buff << s;
            parameters_string += buff.str() + " ";
        }

        string parameters_file_path = out_path;
        string::size_type i = parameters_file_path.rfind('.', parameters_file_path.length());
        if (i != string::npos) {
            parameters_file_path.replace(i + 1, 3, "param");
        }

        std::ofstream out_file(parameters_file_path);
        out_file << parameters_string;
        out_file.close();

        // Create custom mesh resolution file
        string custom_meshing_path = create_custom_meshing_file(resolution);

        // Run command to generate mesh file
        string cmd = "../../isosurface_inflator/isosurface_cli 2D_doubly_periodic " + input_path + " --params '" +
                     parameters_string
                     + "' -m " + custom_meshing_path + " -D inflated.msh -R replicated.msh --cheapPostprocessing " + out_path;

        cout << cmd << endl;

        //system(cmd.c_str());
        execute_command(cmd);
    }

    void generate_topology_and_thickness_info(TReal triangle_side_factor, unsigned num_pillars, double pillar_area_ratio,
                                         TReal thickness_ratio,
                                         Matrix<TReal, 2, Dynamic> &vertices, vector<vector<int>> &edges,
                                         vector<pair<vector<Point>, TReal> > &custom_pairs) {
        double parallelogram_side = 2.0;
        double l = parallelogram_side / 2.0;
        double s = 2 / sqrt(3);
        TReal triangle_side = triangle_side_factor * s * sqrt(3);
        TReal thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars);

        using TPoint = Eigen::Matrix<TReal, 2, 1>;

        unsigned num_upper_intervals = num_pillars - 1;
        TReal upper_spacing = 0.0;
        unsigned num_upper_points = 0;
        if (num_upper_intervals > 0) {
            upper_spacing = (triangle_side - thickness) / (num_pillars - 1);
            num_upper_points = num_pillars / 2 + 1;
        }

        // define important vertices of simplex used to build entire parallelogram structure
        TPoint origin;
        origin << 0, 0;
        TPoint a = origin;
        TPoint b;
        b << l, l * sqrt(3) / 3.0;
        TPoint c;
        c << l, 0;

        // positions can be inferred from the triangle side chosen
        TReal radius = triangle_side * sqrt(3) / 6;
        TPoint aux;
        TPoint f;
        aux << 0, radius;
        f = b - aux;

        // define vertices on simplex based on parameters
        aux(0) = triangle_side / 2.0;
        aux(1) = 0;
        //aux << triangle_side / 2.0, 0;
        TPoint p1 = f - aux;
        TPoint p2 = f;

        Matrix<TReal, 2, Dynamic> simplex_polygon;
        push_back(simplex_polygon, b);
        push_back(simplex_polygon, p1);
        push_back(simplex_polygon, p2);

        vector<Matrix<TReal, 2, Dynamic> > upper_triangles = simplex_polygon_to_whole_parallelogram(simplex_polygon,
                                                                                                     parallelogram_side);

        // Add vertices and edges of the first simplex
        push_back(vertices, b);
        push_back(vertices, p1);
        TPoint p1_offset;
        p1_offset << p1[0] + thickness / 2, p1[1];
        push_back(vertices, p1_offset);

        if ((num_pillars % 2) == 1) {
            push_back(vertices, c);
        }

        // fundamental edges
        edges.push_back({0, 1});
        edges.push_back({1, 2});

        // creating vertices and edges between p1 and p2
        vector<vector<Point> > edges_descriptions;
        vector<Point> vertices_from_p1_to_p2;
        if (num_upper_intervals > 0) {
            create_subdivided_segment_with_constant_spacing(p1_offset, p2, num_upper_points, upper_spacing,
                                                            edges_descriptions);
            vertices_from_p1_to_p2 = extract_vertices_from_edges_descriptions(edges_descriptions);
            add_new_edges(edges_descriptions, vertices, edges);
        } else {
            vertices_from_p1_to_p2.push_back(p2);
            push_back(vertices, p2);
        }

        // adding edge between p2, b and p2,c
        edges.push_back({(int) vertices.cols() - 1, 0});
        if ((num_pillars % 2) == 1) {
            edges.push_back({(int) vertices.cols() - 1, 2});
        }

        vector<Matrix<TReal, 2, Dynamic> > pillar_triangles;
        if (num_upper_intervals >= 0) { // it is always =P
            for (unsigned index = 0; index < vertices_from_p1_to_p2.size(); index++) {
                Point top_point = vertices_from_p1_to_p2[index];
                Point a1, a2, b1, b2;

                if (index == vertices_from_p1_to_p2.size() - 1) {
                    if ((num_pillars % 2) == 0) { // if even number of nodes, last one does not have pillars
                        continue;
                    } else {

                        a1 << top_point[0] - thickness / 2, top_point[1];
                        a2 << top_point[0], top_point[1];

                        b1 << top_point[0] - thickness / 2, 0;
                        b2 << top_point[0], 0;
                    }
                } else {
                    a1 << top_point[0] - thickness / 2, top_point[1];
                    a2 << top_point[0] + thickness / 2, top_point[1];

                    b1 << top_point[0] - thickness / 2, 0;
                    b2 << top_point[0] + thickness / 2, 0;
                }

                Matrix<TReal, 2, 3> triangle_left;
                Matrix<TReal, 2, 3> triangle_right;
                triangle_left << a1, b1, b2;
                triangle_right << a1, b2, a2;

                vector<vector<Point>> new_edges_descriptions = polygon_to_edges_descriptions(triangle_left);
                edges_descriptions.insert(edges_descriptions.end(), new_edges_descriptions.begin(),
                                          new_edges_descriptions.end());
                new_edges_descriptions = polygon_to_edges_descriptions(triangle_right);
                edges_descriptions.insert(edges_descriptions.end(), new_edges_descriptions.begin(),
                                          new_edges_descriptions.end());
                add_new_edges(edges_descriptions, vertices, edges);

                vector<Matrix<TReal, 2, -1, 0, 2, -1>> new_pillar_triangles = simplex_polygon_to_whole_parallelogram(
                        triangle_left, parallelogram_side);
                pillar_triangles.insert(pillar_triangles.end(), new_pillar_triangles.begin(),
                                        new_pillar_triangles.end());
                new_pillar_triangles = simplex_polygon_to_whole_parallelogram(triangle_right, parallelogram_side);
                pillar_triangles.insert(pillar_triangles.end(), new_pillar_triangles.begin(),
                                        new_pillar_triangles.end());
            }
        }

        // create topology in parallelogram
        simplex_to_whole_parallelogram(vertices, edges, parallelogram_side);

        vector<pair<vector<Point>, TReal> > upper_incenters_thickness_pairs = add_polygons_incenters(upper_triangles,
                                                                                                      vertices, edges);
        vector<pair<vector<Point>, TReal> > pillar_triangles_incenters_thickness_pairs = add_polygons_incenters(
                pillar_triangles, vertices, edges);

        custom_pairs = upper_incenters_thickness_pairs;
        custom_pairs.insert(custom_pairs.end(), pillar_triangles_incenters_thickness_pairs.begin(),
                            pillar_triangles_incenters_thickness_pairs.end());
    }

    vector<vector<Point>> create_pillars(vector<Point> segment1, vector<Point> segment2, unsigned n) {
        vector<vector<Point>> edges_descriptions;

        // vectors representing each segment
        Point vector1 = segment1[1] - segment1[0];
        Point vector2 = segment2[1] - segment2[0];

        for (unsigned t=0; t<n; t++) {
            // compute position in segments 1 and 2 that will be linked
            Point pos1 = segment1[0] + vector1 * 1.0 * t / (n - 1);
            Point pos2 = segment2[0] + vector2 * 1.0 * t / (n - 1);

            // create pillars, connecting segments 1 and 2 endpoints
            edges_descriptions.push_back({pos1, pos2});
        }

        return edges_descriptions;
    }

    void create_subdivided_segment(Point segment_start, Point segment_end, unsigned n, vector<vector<Point>> &edges_descriptions) {
        Point vector = segment_end - segment_start;
        Point last_point = segment_start;

        for (unsigned t=1; t<n; t++) {
            Point new_point = segment_start + 1.0 * t / (n - 1) * vector;

            edges_descriptions.push_back({last_point, new_point});
            last_point = new_point;
        }
    }

    void create_subdivided_segment_pillars(Point segment_start, Point segment_end, unsigned n, vector<vector<Point>> &edges_descriptions, TReal thickness) {
        Point vector = segment_end - segment_start;
        Point unit_vector = vector / vector.norm();

        // first point
        Point pillar_start = segment_start - thickness * unit_vector;
        Point pillar_stop = segment_start + thickness * unit_vector;
        edges_descriptions.push_back({pillar_start, segment_start});
        edges_descriptions.push_back({segment_start, pillar_stop});
        Point last_point = pillar_stop;

        for (unsigned t=1; t<n; t++) {
            Point new_point = segment_start + 1.0 * t / (n - 1) * vector;

            pillar_start = new_point - thickness * unit_vector;
            pillar_stop = new_point + thickness * unit_vector;

            edges_descriptions.push_back({last_point, pillar_start});
            edges_descriptions.push_back({pillar_start, new_point});
            edges_descriptions.push_back({new_point, pillar_stop});
            last_point = pillar_stop;
        }

    }

    vector<vector<Point>> create_triangle_edges(Matrix<TReal, 2, 3> triangle, unsigned n, TReal offset=0.0, TReal thickness=0.0) {
        vector<vector<Point>> edges_descriptions;
        Point segment_start, segment_end, segment_start_offset, segment_end_offset;

        for (unsigned i = 0; i<3; i++) {
            segment_start = triangle.col(i % 3);
            segment_end = triangle.col((i + 1) % 3);

            Point unit_vector = (segment_end - segment_start) / (segment_end - segment_start).norm();

            segment_start_offset = segment_start + unit_vector * offset;
            segment_end_offset = segment_end - unit_vector * offset;

            if (thickness > 0.0)
                create_subdivided_segment_pillars(segment_start_offset, segment_end_offset, n, edges_descriptions, thickness);
            else
                create_subdivided_segment(segment_start_offset, segment_end_offset, n, edges_descriptions);

        }

        return edges_descriptions;
    }



    void generate_simpler_topology_and_thickness_info(TReal triangle_side_factor, unsigned num_pillars, double pillar_area_ratio,
                                              TReal thickness_ratio, Matrix<TReal, 2, Dynamic> &vertices, vector<vector<int>> &edges,
                                              vector<pair<vector<Point>, TReal> > &custom_pairs) {
        double parallelogram_side = 2.0;
        double l = parallelogram_side / 2.0;
        double s = 2 / sqrt(3);
        TReal triangle_side = triangle_side_factor * s * sqrt(3);
        TReal thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars) / 2; // using radius as thickness

        using TPoint = Eigen::Matrix<TReal, 2, 1>;

        unsigned num_upper_intervals = num_pillars - 1;
        TReal upper_spacing = 0.0;
        unsigned num_upper_points = 0;
        if (num_upper_intervals > 0) {
            upper_spacing = (triangle_side - thickness) / (num_pillars - 1);
            num_upper_points = num_pillars / 2 + 1;
        }

        // define important vertices of simplex used to build entire parallelogram structure
        TPoint origin;
        origin << 0, 0;
        TPoint triangle_centroid;
        triangle_centroid << l, l * sqrt(3) / 3.0;

        // First side: create pillars of the left side (counter-clock orientation)
        //     /.
        //    /  .
        //   /.....
        TPoint parallelogram_start;
        parallelogram_start << (l + triangle_side / 2.0) / 2.0, (l + triangle_side / 2.0) * sqrt(3) / 2.0;
        TPoint parallelogram_end;
        parallelogram_end << (l - triangle_side / 2.0) / 2.0, (l - triangle_side / 2.0) * sqrt(3) / 2.0;

        TReal radius = triangle_side * sqrt(3) / 6;
        TPoint offset_point;
        offset_point << 0, radius;
        TPoint f = triangle_centroid - offset_point;

        TPoint triangle_start;
        TPoint triangle_end;
        TPoint aux_point;
        aux_point << 0, sqrt(3) * triangle_side / 2.0;
        triangle_start << f + aux_point;
        aux_point << triangle_side / 2.0, 0;
        triangle_end = f - aux_point;

        TPoint unit_vector;
        TPoint unit_vector_triangle;
        unit_vector << (parallelogram_end - parallelogram_start) / (parallelogram_end - parallelogram_start).norm();
        unit_vector_triangle << (triangle_end - triangle_start) / (triangle_end - triangle_start).norm();

        TPoint parallelogram_start_offset, parallelogram_end_offset, triangle_start_offset, triangle_end_offset;
        parallelogram_start_offset = parallelogram_start + thickness * unit_vector;
        parallelogram_end_offset = parallelogram_end - thickness * unit_vector;
        triangle_start_offset = triangle_start + thickness * unit_vector;
        triangle_end_offset = triangle_end - thickness * unit_vector;

        // add new vertices and edges to current sets
        vector<vector<Point>> new_edges_description = create_pillars({parallelogram_start_offset, parallelogram_end_offset},
                                                      {triangle_start_offset, triangle_end_offset}, num_pillars);

        // save pillar nodes
        Matrix<TReal, 2, Dynamic> pillar_nodes;
        for (unsigned index = 0; index < new_edges_description.size(); index++) {
            vector<TPoint> edge_description = new_edges_description[index];

            push_back(pillar_nodes, edge_description[0]);
            push_back(pillar_nodes, edge_description[1]);
        }

        add_new_edges(new_edges_description, vertices, edges);


        // Second side: create pillars of the bottom side (counter-clock orientation)
        //     ..
        //    .  .
        //   ______
        parallelogram_start << l - triangle_side / 2.0, 0;
        parallelogram_end << l + triangle_side / 2.0, 0;

        aux_point << triangle_side / 2.0, 0;
        triangle_start = f - aux_point;
        triangle_end = f + aux_point;

        Matrix<TReal, 2, 2> pillar_example;
        pillar_example << parallelogram_start, triangle_start;

        unit_vector = (parallelogram_end - parallelogram_start) / (parallelogram_end - parallelogram_start).norm();
        unit_vector_triangle = (triangle_end - triangle_start) / (triangle_end - triangle_start).norm();

        parallelogram_start_offset = parallelogram_start + thickness * unit_vector;
        parallelogram_end_offset = parallelogram_end - thickness * unit_vector;
        triangle_start_offset = triangle_start + thickness * unit_vector;
        triangle_end_offset = triangle_end - thickness * unit_vector;

        // add new vertices and edges to current sets
        new_edges_description = create_pillars({parallelogram_start_offset, parallelogram_end_offset},
                                               {triangle_start_offset, triangle_end_offset},
                                               num_pillars);

        // save pillar nodes
        for (unsigned int index = 0; index < new_edges_description.size(); index++) {
            vector<TPoint> edge_description = new_edges_description[index];

            push_back(pillar_nodes, edge_description[0]);
            push_back(pillar_nodes, edge_description[1]);
        }

        add_new_edges(new_edges_description, vertices, edges);


        // Third side: create pillars of the right side (counter-clock orientation)
        //     .\
        //    .  \
        //   .....\

        parallelogram_start << 2 * l - (l - triangle_side / 2.0) / 2, (l - triangle_side / 2.0) * sqrt(3) / 2;
        parallelogram_end << 2 * l - (l + triangle_side / 2.0) / 2, (l + triangle_side / 2.0) * sqrt(3) / 2;

        aux_point << triangle_side / 2.0, 0;
        triangle_start = f + aux_point;
        aux_point << 0, sqrt(3) * triangle_side / 2.0;
        triangle_end = f + aux_point;

        unit_vector = (parallelogram_end - parallelogram_start) / (parallelogram_end - parallelogram_start).norm();
        unit_vector_triangle = (triangle_end - triangle_start) / (triangle_end - triangle_start).norm();

        parallelogram_start_offset = parallelogram_start + thickness * unit_vector;
        parallelogram_end_offset = parallelogram_end - thickness * unit_vector;
        triangle_start_offset = triangle_start + thickness * unit_vector;
        triangle_end_offset = triangle_end - thickness * unit_vector;

        // add new vertices and edges to current sets
        new_edges_description = create_pillars({parallelogram_start_offset, parallelogram_end_offset},
                                                      {triangle_start_offset, triangle_end_offset},
                                                      num_pillars);

        // save nodes on hypotenuse
        Matrix<TReal, 2, Dynamic> hypotenuse_nodes;
        // save pillar nodes
        for (unsigned int index = 0; index < new_edges_description.size(); index++) {
            vector<TPoint> edge_description = new_edges_description[index];

            push_back(hypotenuse_nodes, edge_description[0]);
            push_back(hypotenuse_nodes, edge_description[1]);
        }

        add_new_edges(new_edges_description, vertices, edges);


        // Now, adding edges on the triangle
        aux_point << 0, sqrt(3) * triangle_side / 2.0;
        Point aux_point2;
        aux_point2 << triangle_side / 2.0, 0;
        Matrix<TReal, 2, 3> triangle;
        triangle << f + aux_point, f - aux_point2, f + aux_point2;
        new_edges_description = create_triangle_edges(triangle, num_pillars, thickness, thickness);
        add_new_edges(new_edges_description, vertices, edges);

        // Now, transform to square every vertex we have:
        Matrix<TReal, 2, 2> transformation_matrix;
        transformation_matrix << 1.0, -0.577350269189626,
                0.0, 1.154700538379251;
        Matrix<TReal, 2, Dynamic> resulting_vertices = transformation_matrix * vertices;
        vertices = resulting_vertices;

        // Finally, transpose to origin and reflect through the diagonal 'y = -x'
        offset_point << 1.0, 1.0;
        Matrix<TReal, 2, Dynamic> translated_vertices = vertices.colwise() - offset_point;
        Matrix<TReal, 2, Dynamic> reflected_vertices = reflect(135, translated_vertices);
        vertices = translated_vertices;

        // saving vertices on hypotenuses
        Matrix<TReal, 2, Dynamic> resulting_hypotenuse_nodes = transformation_matrix * hypotenuse_nodes;
        offset_point << 1.0, 1.0;
        Matrix<TReal, 2, Dynamic> translated_hypotenuse_nodes = resulting_hypotenuse_nodes.colwise() - offset_point;
        Matrix<TReal, 2, Dynamic> reflected_hypotenuse_nodes = reflect(135, translated_hypotenuse_nodes);
        Matrix<TReal, Dynamic, Dynamic> complete_hypotenuse_nodes(hypotenuse_nodes.rows(), hypotenuse_nodes.cols() + reflected_hypotenuse_nodes.cols());
        complete_hypotenuse_nodes << translated_hypotenuse_nodes, reflected_hypotenuse_nodes;
        hypotenuse_nodes = complete_hypotenuse_nodes;

        // saving vertices on pillars
        Matrix<TReal, 2, Dynamic> resulting_pillar_nodes = transformation_matrix * pillar_nodes;
        Matrix<TReal, 2, Dynamic> translated_pillar_nodes = resulting_pillar_nodes.colwise() - offset_point;
        Matrix<TReal, 2, Dynamic> reflected_pillar_nodes = reflect(135, translated_pillar_nodes);
        Matrix<TReal, Dynamic, Dynamic> complete_pillar_nodes(pillar_nodes.rows(), pillar_nodes.cols() + reflected_pillar_nodes.cols());
        complete_pillar_nodes << translated_pillar_nodes, reflected_pillar_nodes;
        pillar_nodes = complete_pillar_nodes;


        // transforming pillar example
        pillar_example = transformation_matrix * pillar_example;

        TReal delta_y = abs(pillar_example.col(1)(1) - pillar_example.col(0)(1));
        TReal norm_example = (pillar_example.col(1) - pillar_example.col(0)).norm();
        TReal thickness_correction_factor = delta_y / norm_example;
        cout << "Thickness correction factor: "  << thickness_correction_factor << endl;


        new_edges_description.clear();
        for (unsigned index=0; index < edges.size(); index++) {
            vector<int> edge = edges[index];

            int v1 = edge[0];
            int v2 = edge[1];

            Point reflected_vertex1 = reflected_vertices.col(v1);
            Point reflected_vertex2 = reflected_vertices.col(v2);

            new_edges_description.push_back({reflected_vertex1, reflected_vertex2});
        }

        add_new_edges(new_edges_description, vertices, edges);


        // triangle is also transformed and reflected
        Matrix<TReal, 2, Dynamic> resulting_triangle = transformation_matrix * triangle;
        offset_point << 1.0, 1.0;
        Matrix<TReal, 2, Dynamic> translated_triangle = resulting_triangle.colwise() - offset_point;
        Matrix<TReal, 2, Dynamic> reflected_triangle = reflect(135, translated_triangle);
        Matrix<TReal, 2, Dynamic> triangle_vertices;

        vector<pair<vector<Point>, TReal> > incenter_triangle_pairs = add_polygons_incenters({translated_triangle, reflected_triangle}, vertices, edges);

        // concatenate all custom pairs
        custom_pairs = incenter_triangle_pairs;

        vector<Point> hypotenuse_nodes_vector;
        for (unsigned index=0; index<hypotenuse_nodes.cols(); index++) {
            Point node =  hypotenuse_nodes.col(index);
            hypotenuse_nodes_vector.push_back(node);
        }
        custom_pairs.push_back({hypotenuse_nodes_vector, thickness * sqrt(2)});

        vector<Point> pillar_nodes_vector;
        for (unsigned index=0; index<pillar_nodes.cols(); index++) {
            Point node =  pillar_nodes.col(index);
            pillar_nodes_vector.push_back(node);
        }
        custom_pairs.push_back({pillar_nodes_vector, thickness * thickness_correction_factor});
    }

    TReal get_thickness(TReal min_thickness_ratio, TReal max_thickness_ratio, unsigned index, unsigned num_pillars, TReal total_pillar_region) {
        TReal w_min = min_thickness_ratio * (total_pillar_region / num_pillars);
        TReal w_max = max_thickness_ratio * (total_pillar_region / num_pillars);
        unsigned m = 0;
        TReal angular_coeff, linear_coeff;

        if (num_pillars < 3) {
            return w_min;
        }
        else if (num_pillars % 2 == 1) {
            // odd
            m = (num_pillars + 1) / 2;
            angular_coeff = (w_max - w_min) / (m - 1);
            linear_coeff = w_min;
        }
        else {
            // even
            m = num_pillars / 2;
            angular_coeff = (w_max - w_min) / (m - 1);
            linear_coeff = w_min;
        }

        TReal result;
        if (index > m) {
            result = linear_coeff + (num_pillars - index) * angular_coeff;
        }
        else {
            result = linear_coeff + (index - 1) * angular_coeff;
        }

        cout << "Thickness: " << result << endl;

        return result;
    }

    TReal get_spacing(TReal min_thickness_ratio, TReal max_thickness_ratio, unsigned num_pillars, TReal total_pillar_region) {
        TReal w_min = min_thickness_ratio * (total_pillar_region / num_pillars);
        TReal w_max = max_thickness_ratio * (total_pillar_region / num_pillars);
        int m = 0;
        TReal angular_coeff, linear_coeff;
        TReal sum;

        cout << "w_min: " << w_min << endl;
        cout << "w_max: " << w_max << endl;

        if (num_pillars < 3) {
            TReal spacing = (total_pillar_region - w_min) / (num_pillars - 1) - w_min;
            return spacing;
        }
        else if (num_pillars % 2 == 1) {
            // odd
            m = (num_pillars + 1) / 2;
            angular_coeff = (w_max - w_min) / (m - 1);
            linear_coeff = w_min;

            sum = (2*m - 1) * linear_coeff + angular_coeff * (1 -2*m + pow(m, 2));

            cout << "sum: " << sum << endl;
        }
        else {
            // even
            m = num_pillars / 2;
            angular_coeff = (w_max - w_min) / (m - 1);
            linear_coeff = w_min;

            sum = 2*m * linear_coeff + angular_coeff * (pow(m, 2)-m);
        }

        // computes constant spacing
        TReal result = (total_pillar_region - sum) / (num_pillars - 1);

        cout << "Spacing: " << result << endl;
        return result;
    }

    void generate_auxetic_topology_and_thickness_info(TReal triangle_side_ratio, unsigned num_pillars, TReal pillar_area_ratio,
                                                      TReal min_thickness_ratio, TReal max_thickness_ratio, TReal ninja_factor,
                                                      Matrix<TReal, 2, Dynamic> &vertices,
                                                      vector<vector<int>> &edges, vector<pair<vector<Point>, TReal> > &custom_pairs) {
        vector<vector<Point>> edges_descriptions;
        double parallelogram_side = 3.0;
        double s = parallelogram_side / 3.0;

        TReal p1 = triangle_side_ratio;
        TReal p2 = num_pillars;
        TReal p3 = pillar_area_ratio;
        TReal p4 = min_thickness_ratio;
        TReal p5 = max_thickness_ratio;
        TReal p6 = ninja_factor;

        vector<Matrix<TReal, 2, Dynamic> > triangles, pillar_polygons;

        // define important vertices of simplex used to build entire parallelogram structure
        Point origin, a, b, c;
        origin << 0, 0;
        a << 0, 0;
        b << s/2.0, s * sqrt(3)/2.0;
        c << s, 0;

        // define vertices of triangle
        TReal triangle_side = triangle_side_ratio * s;
        TReal triangle_y_position = s * sqrt(3)/2.0 * (1-triangle_side_ratio);

        Point q1, q2, w, z;
        Point q1_reflected, q2_reflected, w_reflected, z_reflected;
        Point ba_unit, ab_unit;
        q1 << s/2.0 * (1-triangle_side_ratio), triangle_y_position;
        q2 << s/2.0 * (1+triangle_side_ratio), triangle_y_position;

        w  << q2[0] - p1*p3*s, triangle_y_position;
        ba_unit = (b - a) / (b-a).norm();
        z = ba_unit * p1*p3*s*(1-p6) + w;

        q1_reflected << q2[0], -q2[1];
        q2_reflected << q1[0], -q1[1];
        w_reflected  << q2_reflected[0] + p1*p3*s , -triangle_y_position;
        ab_unit = (a-b) / (a-b).norm();
        z_reflected = ab_unit * p1*p3*s*(1-p6) + w_reflected;

        // computing thickness and spacing
        TReal pillar_area = (z - q2).norm();

        // CREATE VERTICES IN THE INITIAL (ORIGINAL) EQUILATERAL TRIANGLE
        create_pillars_with_constant_spacing_and_thickness({q2, z}, {z_reflected, q2_reflected}, min_thickness_ratio, max_thickness_ratio, pillar_area, num_pillars, edges_descriptions, pillar_polygons);
        vector<Matrix<TReal, 2, Dynamic>> original_pillar_polygons;
        std::copy(pillar_polygons.begin(), pillar_polygons.end(), std::back_insert_iterator<vector<Matrix<TReal, 2, Dynamic> > >(original_pillar_polygons));
        edges_descriptions.push_back({z, q1});
        add_new_edges(edges_descriptions, vertices, edges);

        // CREATE BOTTOM-LEFT PART OF THE STRUCTURE
        edges_descriptions.clear();
        rotate_at(b, 60, vertices, edges, edges_descriptions);
        rotate_at(b, 120, vertices, edges, edges_descriptions);

        // do the same with hexagon vertices
        Matrix<TReal, 2, 3> triangle;
        triangle << q1, b, z;
        triangles.push_back(triangle);
        triangle << q2, b, z;
        triangles.push_back(triangle);

        Matrix<TReal, 2, 3> hexagon_vertices_matrix;
        hexagon_vertices_matrix << q1, z, q2;
        vector<vector<int>> hexagon_edges = {{0, 1}, {1, 2}};
        vector<vector<Point>> hexagon_edges_descriptions = {{q1, z}, {z, q2}};
        rotate_at(b, 60, hexagon_vertices_matrix, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(b, 120, hexagon_vertices_matrix, hexagon_edges, hexagon_edges_descriptions);
        vector<Point> hexagon_vertices = extract_vertices_from_edges_descriptions(hexagon_edges_descriptions);
        for (unsigned i = 0; i<6; i++) {
            triangle << b, hexagon_vertices[i], hexagon_vertices[i+1];
            triangles.push_back(triangle);
        }

        // same with pillar polygons
        for (unsigned index = 0; index < original_pillar_polygons.size(); index++) {
            Matrix<TReal, 2, Dynamic> pillar_polygon = original_pillar_polygons[index];

            vector<vector<Point>> empty_edges_descriptions;
            Matrix<TReal, 2, Dynamic> rotated_60 = rotate_at(b, 60, pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_120 = rotate_at(b, 120, pillar_polygon, {}, empty_edges_descriptions);
            pillar_polygons.push_back(rotated_60);
            pillar_polygons.push_back(rotated_120);
        }


        // CREATE TOP-LEFT PART OF THE STRUCTURE
        Point center;
        center <<2*s, s*sqrt(3);
        Matrix<TReal, 2, Dynamic> translated_vertices = translate_at(b, center, vertices, edges, edges_descriptions);
        rotate_at(center, 60, translated_vertices, edges, edges_descriptions);
        rotate_at(center, 120, translated_vertices, edges, edges_descriptions);
        rotate_at(center, 180, translated_vertices, edges, edges_descriptions);
        rotate_at(center, 240, translated_vertices, edges, edges_descriptions);
        rotate_at(center, 300, translated_vertices, edges, edges_descriptions);

        // do the same with hexagon vertices
        hexagon_vertices_matrix << q1, z, q2;
        hexagon_edges = {{0, 1}, {1, 2}};
        hexagon_vertices = {q1, z, q2};
        hexagon_edges_descriptions.clear();
        Matrix<TReal, 2, Dynamic> translated_hexagon_vertices = translate_at(b, center, hexagon_vertices_matrix, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 60, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 120, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 180, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 240, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 300, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 360, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        hexagon_vertices = extract_vertices_from_edges_descriptions(hexagon_edges_descriptions);
        for (unsigned i = 0; i<11; i++) {
            triangle << center, hexagon_vertices[i], hexagon_vertices[i+1];
            triangles.push_back(triangle);
        }
        triangle << center, hexagon_vertices[0], hexagon_vertices[11];
        triangles.push_back(triangle);


        // same with pillar polygons
        for (unsigned index = 0; index < original_pillar_polygons.size(); index++) {
            Matrix<TReal, 2, Dynamic> pillar_polygon = original_pillar_polygons[index];

            vector<vector<Point>> empty_edges_descriptions;

            Matrix<TReal, 2, Dynamic> translated_pillar_polygon = translate_at(b, center, pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_60 = rotate_at(center, 60, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_120 = rotate_at(center, 120, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_180 = rotate_at(center, 180, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_240 = rotate_at(center, 240, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_300 = rotate_at(center, 300, translated_pillar_polygon, {}, empty_edges_descriptions);

            pillar_polygons.push_back(translated_pillar_polygon);
            pillar_polygons.push_back(rotated_60);
            pillar_polygons.push_back(rotated_120);
            pillar_polygons.push_back(rotated_180);
            pillar_polygons.push_back(rotated_240);
            pillar_polygons.push_back(rotated_300);
        }

        // CREATE TOP-RIGHT PART OF THE STRUCTURE
        center << 3*s + s/2, 3.0/2*s*sqrt(3);
        translated_vertices = translate_at(b, center, vertices, edges, edges_descriptions);
        rotate_at(center, 60, translated_vertices, edges, edges_descriptions);
        rotate_at(center, -60, translated_vertices, edges, edges_descriptions);

        // do the same with hexagon vertices
        vector<vector<Point>> empty_edges_descriptions;
        hexagon_vertices_matrix << q1, z, q2;
        hexagon_edges = {{0, 1}, {1, 2}};
        hexagon_vertices = {q1, z, q2};
        hexagon_edges_descriptions.clear();
        translated_hexagon_vertices = translate_at(b, center, hexagon_vertices_matrix, hexagon_edges, empty_edges_descriptions);
        rotate_at(center, -60, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        hexagon_edges_descriptions.push_back({translated_hexagon_vertices.col(0), translated_hexagon_vertices.col(1)});
        hexagon_edges_descriptions.push_back({translated_hexagon_vertices.col(1), translated_hexagon_vertices.col(2)});
        rotate_at(center, 60, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        hexagon_vertices = extract_vertices_from_edges_descriptions(hexagon_edges_descriptions);
        for (unsigned i = 0; i<6; i++) {
            triangle << center, hexagon_vertices[i], hexagon_vertices[i+1];
            triangles.push_back(triangle);
        }

        // same with pillar polygons
        for (unsigned index = 0; index < original_pillar_polygons.size(); index++) {
            Matrix<TReal, 2, Dynamic> pillar_polygon = original_pillar_polygons[index];

            vector<vector<Point>> empty_edges_descriptions;

            Matrix<TReal, 2, Dynamic> translated_pillar_polygon = translate_at(b, center, pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_60 = rotate_at(center, 60, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_300 = rotate_at(center, 300, translated_pillar_polygon, {}, empty_edges_descriptions);

            pillar_polygons.push_back(translated_pillar_polygon);
            pillar_polygons.push_back(rotated_60);
            pillar_polygons.push_back(rotated_300);
        }


        // CREATE BOTTOM PART OF THE STRUCTURE
        center << 2*s, 0.0;
        translated_vertices = translate_at(b, center, vertices, edges, empty_edges_descriptions);
        rotate_at(center, 120, translated_vertices, edges, edges_descriptions);
        rotate_at(center, 180, translated_vertices, edges, edges_descriptions);
        rotate_at(center, 240, translated_vertices, edges, edges_descriptions);

        // do the same with hexagon vertices
        hexagon_vertices = {q1, z, q2};
        hexagon_vertices_matrix << q1, z, q2;
        hexagon_edges = {{0, 1}, {1, 2}};
        hexagon_edges_descriptions.clear();
        translated_hexagon_vertices = translate_at(b, center, hexagon_vertices_matrix, hexagon_edges, empty_edges_descriptions);
        rotate_at(center, 120, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 180, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 240, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        hexagon_vertices = extract_vertices_from_edges_descriptions(hexagon_edges_descriptions);
        for (unsigned i = 0; i<6; i++) {
            triangle << center, hexagon_vertices[i], hexagon_vertices[i+1];
            triangles.push_back(triangle);
        }

        // same with pillar polygons
        for (unsigned index = 0; index < original_pillar_polygons.size(); index++) {
            Matrix<TReal, 2, Dynamic> pillar_polygon = original_pillar_polygons[index];

            vector<vector<Point>> empty_edges_descriptions;

            Matrix<TReal, 2, Dynamic> translated_pillar_polygon = translate_at(b, center, pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_120 = rotate_at(center, 120, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_180 = rotate_at(center, 180, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_240 = rotate_at(center, 240, translated_pillar_polygon, {}, empty_edges_descriptions);

            pillar_polygons.push_back(rotated_120);
            pillar_polygons.push_back(rotated_180);
            pillar_polygons.push_back(rotated_240);
        }


        // CREATE RIGHT PART OF THE STRUCTURE
        center << 3*s + s/2, s*sqrt(3)/2;
        translated_vertices = translate_at(b, center, vertices, edges, empty_edges_descriptions);
        rotate_at(center, 180, translated_vertices, edges, edges_descriptions);
        rotate_at(center, 240, translated_vertices, edges, edges_descriptions);
        rotate_at(center, 300, translated_vertices, edges, edges_descriptions);

        // do the same with hexagon vertices
        hexagon_vertices = {q1, z, q2};
        hexagon_vertices_matrix << q1, z, q2;
        hexagon_edges = {{0, 1}, {1, 2}};
        hexagon_edges_descriptions.clear();
        translated_hexagon_vertices = translate_at(b, center, hexagon_vertices_matrix, hexagon_edges, empty_edges_descriptions);
        rotate_at(center, 180, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 240, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        rotate_at(center, 300, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions);
        hexagon_vertices = extract_vertices_from_edges_descriptions(hexagon_edges_descriptions);
        for (unsigned i = 0; i<6; i++) {
            triangle << center, hexagon_vertices[i], hexagon_vertices[i+1];
            triangles.push_back(triangle);
        }

        // same with pillar polygons
        for (unsigned index = 0; index < original_pillar_polygons.size(); index++) {
            Matrix<TReal, 2, Dynamic> pillar_polygon = original_pillar_polygons[index];

            Matrix<TReal, 2, Dynamic> translated_pillar_polygon = translate_at(b, center, pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_180 = rotate_at(center, 180, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_240 = rotate_at(center, 240, translated_pillar_polygon, {}, empty_edges_descriptions);
            Matrix<TReal, 2, Dynamic> rotated_300 = rotate_at(center, 300, translated_pillar_polygon, {}, empty_edges_descriptions);

            pillar_polygons.push_back(rotated_180);
            pillar_polygons.push_back(rotated_240);
            pillar_polygons.push_back(rotated_300);
        }

        // ADD ALL EDGES DESCRIPTIONS
        add_new_edges(edges_descriptions, vertices, edges);


        // TRANSFORM TO FINAL SQUARE SHAPE in [-1, 1] x [-1, 1] cell
        Matrix<TReal, 2, 2> transformation_matrix;
        transformation_matrix << 1.0, -0.577350269189626,
                0.0, 1.154700538379251;
        transformation_matrix = 2.0/3.0 * transformation_matrix;
        Matrix<TReal, 2, Dynamic> resulting_vertices = transformation_matrix * vertices;
        Point offset;
        offset << 1.0, 1.0;
        vertices = resulting_vertices.colwise() - offset;

        // deal with triangles
        for (unsigned index = 0; index < triangles.size(); index++) {
            Matrix<TReal, 2, Dynamic> triangle = triangles[index];
            Matrix<TReal, 2, Dynamic> resulting_triangle = transformation_matrix * triangle;
            triangle = resulting_triangle.colwise() - offset;
            triangles[index] = triangle;
        }

        vector<pair<vector<Point>, TReal>> incenters_thickness_pairs = add_polygons_incenters(triangles, vertices, edges);

        // deal with pillars
        vector<pair<vector<Point>, TReal>> pillar_nodes_custom_pairs;
        for (unsigned index = 0; index < pillar_polygons.size(); index++) {
            Matrix<TReal, 2, Dynamic> pillar_polygon = pillar_polygons[index];

            Matrix<TReal, 2, Dynamic> resulting_pillar_polygon = transformation_matrix * pillar_polygon;
            pillar_polygon = resulting_pillar_polygon.colwise() - offset;

            TReal nodes_thickness = min_distance_point_line(pillar_polygon.col(0), {pillar_polygon.col(1), pillar_polygon.col(2)});
            vector<Point> pillar_nodes = {pillar_polygon.col(0), pillar_polygon.col(3)};

            pillar_nodes_custom_pairs.push_back({pillar_nodes, nodes_thickness});
        }

        custom_pairs = pillar_nodes_custom_pairs;
        custom_pairs.insert(custom_pairs.end(), incenters_thickness_pairs.begin(), incenters_thickness_pairs.end());
    }



    Matrix<TReal, 2, Dynamic> translate_at(Point origin, Point offset, Matrix<TReal, 2, Dynamic> vertices, vector<vector<int>> edges, vector<vector<Point>> &edges_descriptions) {
        Matrix<TReal, 2, Dynamic> translated_vertices = vertices.colwise() - origin;
        translated_vertices = translated_vertices.colwise() + offset;

        for (unsigned index = 0; index < edges.size(); index++) {
            vector<int> edge = edges[index];

            Point v1 = translated_vertices.col(edge[0]);
            Point v2 = translated_vertices.col(edge[1]);

            edges_descriptions.push_back({v1, v2});
        }

        return translated_vertices;
    }


    Matrix<TReal, 2, Dynamic> rotate_at(Point center_of_rotation, double degrees, Matrix<TReal, 2, Dynamic> vertices, vector<vector<int>> edges, vector<vector<Point>> &edges_descriptions) {
        Matrix<TReal, 2, Dynamic> translated_vertices = vertices.colwise() - center_of_rotation;
        Matrix<TReal, 2, Dynamic> rotated_vertices =
                rotate(degrees, translated_vertices).colwise() + center_of_rotation;

        for (unsigned index = 0; index < edges.size(); index++) {
            vector<int> edge = edges[index];

            Point v1 = rotated_vertices.col(edge[0]);
            Point v2 = rotated_vertices.col(edge[1]);

            edges_descriptions.push_back({v1, v2});
        }

        return rotated_vertices;
    }

    void det_2D(Point a, Point b) {
        return a(0) * b(1) - a(1) * b(0);
    }

    Point edge_intersection(vector<Point> edge1, vector<Point> edge2) {
        Matrix<TReal, 2, 2> M;

        Point x_delta, y_delta;
        x_delta << (edge1[0](0) - edge1[1](0)), (edge2[0](0) - edge2[1](0));
        y_delta << (edge1[0](1) - edge1[1](1)), (edge2[0](1) - edge2[1](1));

        M << x_delta, y_delta;

        TReal div = M.determinant();

        Matrix<TReal, 2, 2> edge1_matrix, edge2_matrix;
        edge1_matrix << edge1[0], edge1[1];
        edge2_matrix << edge2[0], edge2[1];

        Point d;
        d << edge1_matrix.determinant(), edge2_matrix.determinant();

        Matrix<TReal, 2, 2> x_matrix, y_matrix;
        x_matrix << d, x_delta;
        y_matrix << d, y_delta;
        TReal x = x_matrix.determinant() / div;
        TReal y = y_matrix.determinant() / div;

        Point result;
        result << x, y;
        return result;
    }

    vector<Point> create_pillar_edge(vector<Point> edge) {
        Point left, right;
        left << -1.0, 0.0;
        right << 1.0, 0.0;

        vector<Point> crossing_edge = {left, right};
        Point edge_intersection_point = edge_intersection(edge, crossing_edge);
        vector<Point> resulting_edge = {edge[0], edge_intersection_point};

        return resulting_edge;
    }

    double get_pillar_area(TReal triangle_side_ratio, unsigned num_pillars, TReal pillar_area_ratio,
                           TReal min_thickness_ratio, TReal max_thickness_ratio, TReal ninja_factor) {
        double parallelogram_side = 3.0;
        double s = parallelogram_side / 3.0;
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

        w  << q2[0] - triangle_side_ratio*pillar_area_ratio*s, triangle_y_position;
        ba_unit = (b - a) / (b-a).norm();
        z = ba_unit * triangle_side_ratio*pillar_area_ratio*s*(1-ninja_factor) + w;

        q1_reflected << q2[0], -q2[1];
        q2_reflected << q1[0], -q1[1];
        w_reflected  << q2_reflected[0] + triangle_side_ratio*pillar_area_ratio*s , -triangle_y_position;
        ab_unit = (a-b) / (a-b).norm();
        z_reflected = ab_unit * triangle_side_ratio*pillar_area_ratio*s*(1-ninja_factor) + w_reflected;

        TReal pillar_area = (z - q2).norm();

        return pillar_area;
    }

    void create_pillars_with_constant_spacing_and_thickness(vector<Point> line1, vector<Point> line2, TReal min_thickness_ratio, TReal max_thickness_ratio, TReal pillar_area, unsigned num_pillars, vector<vector<Point>> &edges_descriptions, vector<Matrix<TReal, 2, Dynamic>> &pillar_polygons) {
        Point unit_vector = (line1[1] - line1[0]) / (line1[1] - line1[0]).norm();

        // initial points for each side
        Point next_point1 = line1[0];
        Point next_point2 = line2[0];

        Point end_pillar1;

        for (unsigned index = 0; index < num_pillars; index++) {
            TReal thickness = get_thickness(min_thickness_ratio,max_thickness_ratio, index+1, num_pillars, pillar_area);
            TReal gap = get_spacing(min_thickness_ratio, max_thickness_ratio, num_pillars, pillar_area);

            Point beginning_pillar1 = next_point1;
            Point midpoint_pillar1 = beginning_pillar1 + thickness * unit_vector / 2;
            end_pillar1 = beginning_pillar1 + thickness * unit_vector;

            Point beginning_pillar2 = next_point2;
            Point midpoint_pillar2 = beginning_pillar2 + thickness * unit_vector / 2;
            Point end_pillar2 = beginning_pillar2 + thickness * unit_vector;

            vector<Point> top_edge1 = {beginning_pillar1, midpoint_pillar1};
            vector<Point> top_edge2 = {midpoint_pillar1, end_pillar1};

            vector<Point> pillar_edge = {midpoint_pillar1, midpoint_pillar2};
            vector<Point> resulting_pillar_edge = create_pillar_edge(pillar_edge);

            vector<vector<Point>> pillar_edges = {top_edge1, top_edge2, resulting_pillar_edge};

            Matrix<TReal, 2, 4> polygon;
            polygon << midpoint_pillar1, beginning_pillar1, beginning_pillar2, resulting_pillar_edge[1];
            pillar_polygons.push_back(polygon);

            next_point1 = end_pillar1 + gap * unit_vector;
            next_point2 = end_pillar2 + gap * unit_vector;

            if (index < num_pillars - 1) {
                vector<Point> gap_edge = {end_pillar1, next_point1};
                pillar_edges.push_back(gap_edge);
            }

            for (unsigned j = 0; j < pillar_edges.size(); j++) {
                pillar_edge = pillar_edges[j];
                if (!pillar_edge.empty())
                    edges_descriptions.push_back(pillar_edge);
            }
        }

        // need to add final edge connecting next_point1 to end of line1 (only if they differ)
        //vector<Point> final_edge = {end_pillar1, line1[1]};
        //edges_descriptions.push_back(final_edge);

    }
};
#endif //MICROSTRUCTURES_HEXLIB_H
