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

typedef Vector2d Point;

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


    for (unsigned index=0; index<vertices.cols(); index++) {
        Point vertex = vertices.col(index);
        out_file << "v " << std::fixed << std::setprecision(17) << vertex(0) << " " << vertex(1) << " 0" << endl;
    }

    for (unsigned index=0; index<edges.size(); index++) {
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
    unsigned grid_size = pow(2,coarsening) * resolution;
    meshing_opts["marchingSquaresGridSize"] = grid_size;
    meshing_opts["marchingSquaresCoarsening"] = coarsening;

    // write to file
    std::ofstream out_file(custom_name);
    out_file << meshing_opts << std::endl;

    return custom_name;
}

Matrix<double, 2, Dynamic> reflect(double angle, Matrix<double, 2, Dynamic> vertices) {
    double theta = M_PI * (angle / 180);
    Matrix<double, 2, 2> Rl;
    Rl << cos(2 * theta), sin(2 * theta),
            sin(2 * theta), -cos(2 * theta);

    Matrix<double, 2, Dynamic> reflected_vertices = Rl * vertices;

    return reflected_vertices;
}

Matrix<double, 2, Dynamic> rotate(double angle, Matrix<double, 2, Dynamic> vertices) {
    double theta = M_PI * (angle / 180);
    Matrix<double, 2, 2> R;
    R << cos(theta), -sin(theta),
            sin(theta), cos(theta);

    Matrix<double, 2, Dynamic> rotated_vertices = R * vertices;

    return rotated_vertices;
}

Matrix<double, 2, Dynamic> push_back(Matrix<double, 2, Dynamic> &matrix, Vector2d column) {
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
        if (equal_vertices(v, vertex))
        {
            vertex_position = index;
            break;
        }
    }
    return vertex_position;
}

int find_vertex_in_list(Matrix<double, 2, Dynamic> &vertices, Point v) {
    int vertex_position = -1;

    for (unsigned int index = 0; index < vertices.cols(); index++) {
        Point vertex = vertices.col(index);
        if (equal_vertices(v, vertex))
        {
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

    for (unsigned index=0; index<edges.size(); index++) {
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


double min_distance_point_line(Point point, vector<Point> line) {
    double x0 = point[0];
    double y0 = point[1];
    double x1 = line[0](0);
    double y1 = line[0](1);
    double x2 = line[1](0);
    double y2 = line[1](1);

    double line_norm = sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));

    double distance = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / line_norm;

    return distance;
}


Point triangle_incenter(Matrix<double, 2, Dynamic> triangle) {
    Point p_a = triangle.col(0);
    Point p_b = triangle.col(1);
    Point p_c = triangle.col(2);

    double a = (p_b - p_c).norm();
    double b = (p_a - p_c).norm();
    double c = (p_a - p_b).norm();

    Point incenter = (a * p_a + b * p_b + c * p_c) / (a + b + c);

    return incenter;
}

vector<vector<Point>> polygon_to_edges_descriptions(Matrix<double, 2, -1, 0, 2, -1> polygon) {
    vector<vector<Point>> edges_descriptions;

    for (unsigned i=0; i<polygon.cols(); i++) {
        Point p1 = polygon.col(i);
        Point p2 = polygon.col((i + 1) % polygon.cols());

        edges_descriptions.push_back({p1, p2});
    }

    return edges_descriptions;
}

void add_new_edges(vector< vector<Point> > &edges_descriptions, Matrix<double, 2, Dynamic> &current_vertices, vector<vector<int>> &current_edges) {

    for (unsigned index=0; index < edges_descriptions.size(); index++) {
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

vector< pair<vector<Point>, double> > add_polygons_incenters(vector<Matrix<double, 2, Dynamic> > polygons, Matrix<double, 2, Dynamic> &vertices, vector<vector<int>> &edges) {
    vector< pair<vector<Point>, double> > incenters_thickness_pairs;

    for (unsigned index=0; index<polygons.size(); index++) {
        Matrix<double, 2, Dynamic> pol = polygons[index];

        Point incenter = triangle_incenter(pol);
        double thickness = min(min(min_distance_point_line(incenter, {pol.col(0), pol.col(1)}),
                                   min_distance_point_line(incenter, {pol.col(1), pol.col(2)})),
                               min_distance_point_line(incenter, {pol.col(2), pol.col(0)}));

        incenters_thickness_pairs.push_back({{incenter}, thickness});

        vector<vector<Point>> edges_descriptions;

        for (int i=0; i<pol.cols(); i++) {
            Point a = pol.col(i);
            Point b = incenter;
            edges_descriptions.push_back({a, b});
        }

        add_new_edges(edges_descriptions, vertices, edges);
    }

    return incenters_thickness_pairs;
}

vector<Point> extract_vertices_from_edges_descriptions(vector< vector<Point> >edges_descriptions) {
    vector<Point> extracted_vertices;

    for (unsigned int index=0; index<edges_descriptions.size(); index++)
    {
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


void create_subdivided_segment_with_constant_spacing(Point segment_start, Point segment_end, int n, double spacing, vector< vector<Point> > &edges_descriptions) {
    Vector2d vector = (segment_end - segment_start);
    Vector2d direction = vector / vector.norm();

    Point last_point = segment_start;

    for (int t=1; t<(n-1); t++) {
        Point new_point = segment_start + spacing * t * direction;
        edges_descriptions.push_back({last_point, new_point});
        last_point = new_point;
    }

    edges_descriptions.push_back({last_point, segment_end});
}

void add_new_vertices_and_edges(Matrix<double, 2, Dynamic> reflected_vertices, vector<vector<int>> reflected_edges, Matrix<double, 2, Dynamic> &vertices, vector<vector<int>> &edges){
    vector< vector<Point> > new_edges_description;

    for (vector<vector<int>>::iterator it = reflected_edges.begin(); it!= reflected_edges.end(); it++) {
        vector<int> edge = *it;

        int v1 = edge[0];
        int v2 = edge[1];

        Point reflected_vertex1 = reflected_vertices.col(v1);
        Point reflected_vertex2 = reflected_vertices.col(v2);

        new_edges_description.push_back({reflected_vertex1, reflected_vertex2});
    }

    add_new_edges(new_edges_description, vertices, edges);
}

void simplex_to_whole_parallelogram(Matrix<double, 2, Dynamic> &vertices, vector<vector<int>> &edges, double parallelogram_side){
    double l = parallelogram_side / 2;

    // Operation 1: reflect against 30 degrees line
    Matrix<double, 2, Dynamic> reflected_vertices = reflect(30, vertices);
    vector<vector<int>> reflected_edges;
    std::copy(edges.begin(), edges.end(), std::back_insert_iterator< vector< vector<int> > >(reflected_edges));
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges);
    Matrix<double, 2, Dynamic> vertices_step_1;
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
    std::copy(edges.begin(), edges.end(), std::back_insert_iterator< vector< vector<int> > >(reflected_edges));
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges);
    vertices = vertices.colwise() + offset_point;

    // Operation 3: rotate 240 degrees and plug on top of current structure
    offset_point << l, sqrt(3) / 3;
    Matrix<double, 2, Dynamic> translated_vertices = vertices_step_1.colwise() - offset_point;
    Matrix<double, 2, Dynamic> rotated_vertices = rotate(240, translated_vertices);
    vector<vector<int>> rotated_edges;
    std::copy(edges_step_1.begin(), edges_step_1.end(), std::back_insert_iterator< vector< vector<int> > >(rotated_edges));
    rotated_vertices = rotated_vertices.colwise() + offset_point;
    add_new_vertices_and_edges(rotated_vertices, rotated_edges, vertices, edges);

    // Operation 4:
    offset_point << 2 * l, 0;
    translated_vertices = vertices.colwise() - offset_point;
    reflected_vertices = reflect(120, translated_vertices);
    reflected_edges.clear();
    std::copy(edges.begin(), edges.end(), std::back_insert_iterator< vector< vector<int> > >(reflected_edges));
    reflected_vertices = reflected_vertices.colwise() + offset_point;
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges);

    // Operation 5: transform parallelogram to square
    Matrix<double, 2, 2> transformation_matrix;
    transformation_matrix << 1.0, -0.577350269189626,
            0.0, 1.154700538379251;
    Matrix<double, 2, Dynamic> resulting_vertices = transformation_matrix * vertices;
    vertices = resulting_vertices;

    // Finally, transpose to origin
    offset_point << l, l;
    vertices = vertices.colwise() - offset_point;
}


vector<Matrix<double, 2, Dynamic> > simplex_polygon_to_whole_parallelogram(Matrix<double, 2, Dynamic> simplex_polygon, double parallelogram_side) {
    double l = parallelogram_side / 2;
    Matrix<double, 2, Dynamic> vertices = simplex_polygon;
    vector<Matrix<double, 2, Dynamic> > polygons;
    polygons.emplace_back(simplex_polygon);

    // Operation 1: reflect against 30 degrees line
    Matrix<double, 2, Dynamic> reflected_vertices = reflect(30, vertices);
    polygons.emplace_back(reflected_vertices);  // polygons now has 2 pols
    vector<Matrix<double, 2, Dynamic> > polygons_step_1;
    std::copy(polygons.begin(), polygons.end(),
              std::back_insert_iterator<vector<Matrix<double, 2, Dynamic> > >(polygons_step_1));

    // Operation 2: reflect against y axis at x = l
    vector<Matrix<double, 2, Dynamic> > polygons_copy;
    std::copy(polygons.begin(), polygons.end(), std::back_insert_iterator<vector<Matrix<double, 2, Dynamic> > >(polygons_copy));
    for (std::vector<Matrix<double, 2, Dynamic> >::iterator it = polygons_copy.begin(); it != polygons_copy.end(); ++it) {
        Point offset_point;
        offset_point << l, 0;
        Matrix<double, 2, Dynamic> pol_vertices = (*it).colwise() - offset_point;
        reflected_vertices = reflect(90, pol_vertices);
        reflected_vertices = reflected_vertices.colwise() + offset_point;
        polygons.emplace_back(reflected_vertices); // polygons now has 4 pols
    }

    // Operation 3: rotate 240 degrees and plug on top of current structure
    polygons_copy.clear();
    std::copy(polygons_step_1.begin(), polygons_step_1.end(),
              std::back_insert_iterator<vector<Matrix<double, 2, Dynamic> > >(polygons_copy));
    for (std::vector<Matrix<double, 2, Dynamic> >::iterator it = polygons_copy.begin(); it != polygons_copy.end(); ++it) {
        Point offset_point;
        offset_point << l, sqrt(3) / 3;
        Matrix<double, 2, Dynamic> translated_vertices = (*it).colwise() - offset_point;
        Matrix<double, 2, Dynamic> rotated_vertices = rotate(240, translated_vertices);
        rotated_vertices = rotated_vertices.colwise() + offset_point;
        polygons.emplace_back(rotated_vertices); // polygons now has 6 pols
    }

    // Operation 4:
    polygons_copy.clear();
    std::copy(polygons.begin(), polygons.end(), std::back_insert_iterator<vector<Matrix<double, 2, Dynamic> > >(polygons_copy));
    for (std::vector<Matrix<double, 2, Dynamic> >::iterator it = polygons_copy.begin(); it != polygons_copy.end(); ++it) {
        Point offset_point;
        offset_point << 2 * l, 0;
        Matrix<double, 2, Dynamic> translated_vertices = (*it).colwise() - offset_point;
        reflected_vertices = reflect(120, translated_vertices);
        reflected_vertices = reflected_vertices.colwise() + offset_point;
        polygons.emplace_back(reflected_vertices); // polygons now has 12 pols
    }

    // Operation 5: transform parallelogram to square
    Matrix<double, 2, 2> transformation_matrix;
    transformation_matrix << 1.0, -0.577350269189626,
            0.0, 1.154700538379251;
    for (unsigned int index=0; index< polygons.size(); index++) {
        Matrix<double, 2, Dynamic> vertices = polygons[index];
        Matrix<double, 2, Dynamic> resulting_vertices = transformation_matrix * vertices;
        Point offset_point;
        offset_point << l, l;
        polygons[index] = resulting_vertices.colwise() - offset_point;
    }

    return polygons;
}

bool is_zero(double number) {
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

    for (unsigned index=0; index<vertices.size(); index++) {
        Point vertex = vertices[index];
        Point independent_vertex = independent_vertex_position(vertex);

        if (equal_vertices(vertex, independent_vertex))
            independent_vertices.push_back(vertex);
    }

    return independent_vertices;
}

vector<double> generate_position_parameters(vector<Point> vertices) {
    vector<double> parameters;

    for (unsigned i=0; i<vertices.size(); i++) {
        Point vertex = vertices[i];

        // first, normalize vertex
        Point normalized_vertex;
        normalized_vertex(0) = (vertex(0) + 1.0) / 2;
        normalized_vertex(1) = (vertex(1) + 1.0) / 2;

        // then, verify if each coordinate is free and, if so, compute parameter
        if (!is_zero(normalized_vertex(0)) &&  !is_zero(normalized_vertex(0) - 1.0)) {
            parameters.push_back(normalized_vertex(0));
        }

        if (!is_zero(normalized_vertex(1)) &&  !is_zero(normalized_vertex(1) - 1.0)) {
            parameters.push_back(normalized_vertex(1));
        }

        // cause z is always at 0 coordinate
        parameters.push_back(0.5);
    }

    return parameters;
}

vector<double> generate_parameters(vector<Point> independent_vertices, double default_thickness, double default_blending, vector<pair<vector<Point>, double>> custom_thickness_pairs) {
    vector<double> parameters = generate_position_parameters(independent_vertices);

    for (unsigned i=0; i<independent_vertices.size(); i++) {
        parameters.push_back(default_thickness);
    }

    for (unsigned i=0; i<independent_vertices.size(); i++) {
        parameters.push_back(default_blending);
    }

    size_t offset_thickeness = parameters.size() - 2 * independent_vertices.size();

    // Now, for the extra parameters, find vertices and apply customized thickness
    for (unsigned index=0; index<custom_thickness_pairs.size(); index++) {
        pair<vector<Point>, double> thickness_pair = custom_thickness_pairs[index];

        vector<Point> custom_nodes = thickness_pair.first;
        double custom_thickness = thickness_pair.second;

        for (vector<Point>::iterator node_it = custom_nodes.begin(); node_it != custom_nodes.end(); node_it++) {
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

void inflate_hexagonal_box(string input_path, double default_thickness, double default_blending, string out_path, vector<pair<vector<Point>, double>> custom_thickness_pairs = vector<pair<vector<Point>, double>>(), unsigned resolution=64) {

    // Discover the vertices
    vector<Point> vertices;


    vector<Point> independent_vertices = extract_independent_vertices(vertices);

    vector<double> parameters = generate_parameters(independent_vertices, default_thickness, default_blending, custom_thickness_pairs);

    // Write parameters to file
    string parameters_string;
    for (float s : parameters) {
        std::ostringstream buff;
        buff<<s;
        parameters_string += buff.str() + " ";
    }

    string parameters_file_path = out_path;
    string::size_type i = parameters_file_path.rfind('.', parameters_file_path.length());
    if (i != string::npos) {
        parameters_file_path.replace(i+1, 3, "param");
    }

    std::ofstream out_file(parameters_file_path);
    out_file << parameters_string;
    out_file.close();

    // Create custom mesh resolution file
    string custom_meshing_path = create_custom_meshing_file(resolution);

    // Run command to generate mesh file
    string cmd = "../../isosurface_inflator/isosurface_cli 2D_doubly_periodic " + input_path + " --params '" + parameters_string
          + "' -m " + custom_meshing_path + " -D inflated.msh -R replicated.msh " + out_path;

    cout << cmd << endl;

    system(cmd.c_str());
    execute_command(cmd);
}

void generate_topology_and_thickness_info(double triangle_side_factor, unsigned num_pillars, double pillar_area_ratio, double thickness_ratio,
                                          Matrix<double, 2, Dynamic> &vertices, vector<vector<int>> &edges, vector<pair<vector<Point>, double> > &custom_pairs) {
    double parallelogram_side = 2.0;
    double l = parallelogram_side / 2.0;
    double s = 2 / sqrt(3);
    double triangle_side = triangle_side_factor * s * sqrt(3);
    double thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars);

    unsigned num_upper_intervals = num_pillars - 1;
    double upper_spacing = 0.0;
    unsigned num_upper_points = 0;
    if (num_upper_intervals > 0) {
        upper_spacing = (triangle_side - thickness) / (num_pillars - 1);
        num_upper_points = num_pillars / 2 + 1;
    }

    // define important vertices of simplex used to build entire parallelogram structure
    Point origin;
    origin << 0, 0;
    Point a = origin;
    Point b;
    b << l, l * sqrt(3) / 3.0;
    Point c;
    c << l, 0;

    // positions can be inferred from the triangle side chosen
    double radius = triangle_side * sqrt(3) / 6;
    Point aux;
    Point f;
    aux << 0, radius;
    f = b - aux;

    // define vertices on simplex based on parameters
    aux << triangle_side / 2.0, 0;
    Point p1 = f - aux;
    Point p2 = f;

    Matrix<double, 2, Dynamic> simplex_polygon;
    push_back(simplex_polygon, b);
    push_back(simplex_polygon, p1);
    push_back(simplex_polygon, p2);

    vector<Matrix<double, 2, Dynamic> > upper_triangles = simplex_polygon_to_whole_parallelogram(simplex_polygon,
                                                                                                 parallelogram_side);

    // Add vertices and edges of the first simplex
    push_back(vertices, b);
    push_back(vertices, p1);
    Point p1_offset;
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

    vector<Matrix<double, 2, Dynamic> > pillar_triangles;
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

            Matrix<double, 2, 3> triangle_left;
            Matrix<double, 2, 3> triangle_right;
            triangle_left << a1, b1, b2;
            triangle_right << a1, b2, a2;

            vector<vector<Point>> new_edges_descriptions = polygon_to_edges_descriptions(triangle_left);
            edges_descriptions.insert(edges_descriptions.end(), new_edges_descriptions.begin(),
                                      new_edges_descriptions.end());
            new_edges_descriptions = polygon_to_edges_descriptions(triangle_right);
            edges_descriptions.insert(edges_descriptions.end(), new_edges_descriptions.begin(),
                                      new_edges_descriptions.end());
            add_new_edges(edges_descriptions, vertices, edges);

            vector<Matrix<double, 2, -1, 0, 2, -1>> new_pillar_triangles = simplex_polygon_to_whole_parallelogram(
                    triangle_left, parallelogram_side);
            pillar_triangles.insert(pillar_triangles.end(), new_pillar_triangles.begin(), new_pillar_triangles.end());
            new_pillar_triangles = simplex_polygon_to_whole_parallelogram(triangle_right, parallelogram_side);
            pillar_triangles.insert(pillar_triangles.end(), new_pillar_triangles.begin(), new_pillar_triangles.end());
        }
    }

    // create topology in parallelogram
    simplex_to_whole_parallelogram(vertices, edges, parallelogram_side);

    vector<pair<vector<Point>, double> > upper_incenters_thickness_pairs = add_polygons_incenters(upper_triangles,
                                                                                                  vertices, edges);
    vector<pair<vector<Point>, double> > pillar_triangles_incenters_thickness_pairs = add_polygons_incenters(
            pillar_triangles, vertices, edges);

    custom_pairs = upper_incenters_thickness_pairs;
    custom_pairs.insert(custom_pairs.end(), pillar_triangles_incenters_thickness_pairs.begin(),
                        pillar_triangles_incenters_thickness_pairs.end());
}

#endif //MICROSTRUCTURES_HEXLIB_H
