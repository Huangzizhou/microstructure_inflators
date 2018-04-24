import math
import os
import re
import json
import uuid
from subprocess import call, check_output

import numpy as np

Tolerance = 1e-10

script_directory = os.getcwd()

def create_wire(vertices, edges, out_wire):
    out_file = open(out_wire, 'w')
    for vertex in vertices:
        out_file.write("v " + str(vertex[0]) + " " + str(vertex[1]) + " 0\n")
    for edge in edges:
        out_file.write("l " + str(edge[0] + 1) + " " + str(edge[1] + 1) + "\n")
    out_file.close()


def create_diamond(initial, end, thickness):
    v = np.array(end) - initial
    v_rot = rotate(90, np.array([v]))
    u_rot = (v_rot / np.linalg.norm(v_rot))[0]

    midpoint = (np.array(end) + initial) / 2.0

    left = midpoint - u_rot * thickness / 2.0
    right = midpoint + u_rot * thickness / 2.0

    return [initial, left, end, right]


def min_distance_to_other_vertices(vertices, centroid_index):
    min_distance = 2.0

    for index, vertex in enumerate(vertices):
        if centroid_index == index:
            continue

        distance = np.linalg.norm(np.array(vertex) - np.array(vertices[centroid_index]))

        if distance < min_distance:
            min_distance = distance

    return min_distance


def create_triangle_edges(triangle, n, offset=0.0, thickness=0.0):
    edges_descriptions = []

    for i in range(0, 3):
        segment_start = triangle[i % 3]
        segment_end = triangle[(i + 1) % 3]

        unit_vector = (segment_end - segment_start) / np.linalg.norm(segment_end - segment_start)

        if n == 1:
            middle_point = segment_start + (segment_end - segment_start) / 2
            segment_start_offset = middle_point - thickness * unit_vector
            segment_end_offset = middle_point + thickness * unit_vector

            edges_descriptions.append([segment_start, segment_start_offset])
            edges_descriptions.append([segment_start_offset, middle_point])
            edges_descriptions.append([middle_point, segment_end_offset])
            edges_descriptions.append([segment_end_offset, segment_end])

            continue

        segment_start_offset = segment_start + unit_vector * offset
        segment_end_offset = segment_end - unit_vector * offset

        if thickness > 0.0:
            create_subdivided_segment_pillars(segment_start_offset, segment_end_offset, n, edges_descriptions, thickness)
        else:
            create_subdivided_segment(segment_start_offset, segment_end_offset, n, edges_descriptions)

    return edges_descriptions


def create_subdivided_segment(segment_start, segment_end, n, edges_descriptions):
    vector = segment_end - segment_start
    last_point = segment_start
    for t in range(1, n):
        new_point = segment_start + 1.0 * t / (n - 1) * vector

        edges_descriptions.append([last_point, new_point])
        last_point = new_point


def create_subdivided_segment_pillars(segment_start, segment_end, n, edges_descriptions, thickness):
    vector = segment_end - segment_start
    unit_vector = vector / np.linalg.norm(vector)

    # first point
    pillar_start = segment_start - thickness * unit_vector
    pillar_stop = segment_start + thickness * unit_vector
    edges_descriptions.append([pillar_start, segment_start])
    edges_descriptions.append([segment_start, pillar_stop])
    last_point = pillar_stop

    for t in range(1, n):
        new_point = segment_start + 1.0 * t / (n - 1) * vector

        pillar_start = new_point - thickness * unit_vector
        pillar_stop = new_point + thickness * unit_vector

        edges_descriptions.append([last_point, pillar_start])
        edges_descriptions.append([pillar_start, new_point])
        edges_descriptions.append([new_point, pillar_stop])
        last_point = pillar_stop


def create_subdivided_segment_with_constant_spacing(segment_start, segment_end, n, spacing, edges_descriptions):
    vector = (segment_end - segment_start)
    direction = vector / np.linalg.norm(vector)
    last_point = segment_start
    for t in range(1, n - 1):
        new_point = segment_start + spacing * t * direction

        edges_descriptions.append([last_point, new_point])
        last_point = new_point

    edges_descriptions.append([last_point, segment_end])


def create_pillars(segment1, segment2, n):
    edges_descriptions = []

    # vectors representing each segment
    vector1 = segment1[1] - segment1[0]
    vector2 = segment2[1] - segment2[0]

    if n == 1:
        pos1 = segment1[0] + vector1 / 2
        pos2 = segment2[0] + vector2 / 2
        edges_descriptions.append([pos1, pos2])

        return edges_descriptions

    for t in range(0, n):
        # compute position in segments 1 and 2 that will be linked
        pos1 = segment1[0] + vector1 * 1.0 * t / (n - 1)
        pos2 = segment2[0] + vector2 * 1.0 * t / (n - 1)

        # create pillars, connecting segments 1 and 2 endpoints
        edges_descriptions.append([pos1, pos2])

    return edges_descriptions


def add_new_edges(edges_descriptions, current_vertices, current_edges):
    for index, edge_description in enumerate(edges_descriptions):
        vertex1 = edge_description[0]
        vertex2 = edge_description[1]

        position1 = find_vertex_in_list(current_vertices, vertex1)
        if position1 == -1:
            position1 = len(current_vertices)
            current_vertices.append(vertex1)

        position2 = find_vertex_in_list(current_vertices, vertex2)
        if position2 == -1:
            position2 = len(current_vertices)
            current_vertices.append(vertex2)

        new_edge = [position1, position2]
        if position1 == position2:
            print("Warning: trying to add edge where both vertices are the same")
            continue

        edge_position = find_edge_in_list(current_edges, new_edge)
        if edge_position == -1:
            current_edges.append(new_edge)


def find_edge_in_list(edges, e):
    new_v1 = e[0]
    new_v2 = e[1]

    edge_position = -1

    for index, edge in enumerate(edges):
        v1 = edge[0]
        v2 = edge[1]

        if (new_v1 == v1) and (new_v2 == v2):
            edge_position = index
            break

        if (new_v1 == v2) and (new_v2 == v1):
            edge_position = index
            break

    return edge_position


def equal_vertices(v1, v2):
    if (abs(v1[0] - v2[0]) < Tolerance) and (abs(v1[1] - v2[1]) < Tolerance):
        return True

    return False


def find_vertex_in_list(vertices, v):
    vertex_position = -1

    for index, vertex in enumerate(vertices):

        if equal_vertices(v, vertex):
            vertex_position = index
            break

    return vertex_position


def inflate_hexagonal_box(input_path, vertices_thickeness, out_path, vertices=[], triangle_centroids=[],
                          hypotenuse_nodes=[]):
    # discover number of vertices
    num_vertices = 0
    with open(input_path, 'r') as input_wire:
        wire_content = input_wire.readlines()
        for line in wire_content:
            if line.startswith('v '):
                num_vertices += 1

    # run isosurface_inflator first time to obtain default parameters
    cwd = os.getcwd()
    with open('default-parameters.txt', 'w') as out_log:
        cmd = [cwd + '/../../isosurface_inflator/isosurface_cli', '2D_doubly_periodic', input_path]
        call(cmd, stdout=out_log)

    with open('default-parameters.txt', 'r') as out_log:
        # read outlog to obtain default parameters and modify them
        content = out_log.readlines()
        parameters = []
        for line_index in range(1, len(content)):
            for element in content[line_index].split():
                # print element
                parameters.append(float(element))

    offset_thickeness = len(parameters) - 2 * num_vertices
    for i in range(offset_thickeness, offset_thickeness + num_vertices):
        parameters[i] = vertices_thickeness

    # find vertices on triangle hypotenuse and adjust their thickness
    for node in hypotenuse_nodes:
        index = find_vertex_in_list(vertices, node)
        if index == -1:
            print("Warning: no vertex in this position", node)
            continue
        else:
            # has to adjust the size because original equilateral triangle in the parallelogram cells is not equilateral
            # after transformed to square
            node_thickness = float(vertices_thickeness) * math.sqrt(2)
            parameters[offset_thickeness + index] = node_thickness

    # find vertices on triangle hypotenuse and adjust their thickness
    for centroid in triangle_centroids:
        index = find_vertex_in_list(vertices, centroid)
        min_distance = min_distance_to_other_vertices(vertices, index)
        parameters[offset_thickeness + index] = max(float(vertices_thickeness), 1.0 * min_distance)

    parameters_string = ' '.join(str(param) for param in parameters)

    cmd = [cwd + '/../../isosurface_inflator/isosurface_cli', '2D_doubly_periodic', input_path, '--params',
           parameters_string, '-m', 'refined-meshing_opts.json', '-D', 'inflated.msh', '-R', 'replicated.msh', out_path]
    # print cmd
    call(cmd)


def is_zero(number):
    if (abs(number) < Tolerance):
        return True

    return False


def independent_vertex_position(vertex):
    v = list(vertex)

    if is_zero(abs(v[0] - 1.0)):
        v[0] = -1.0;

    if is_zero(abs(v[1] - 1.0)):
        v[1] = -1.0;

    return v


def extract_independent_vertices(vertices):
    num_indep_vertices = len(vertices)
    num_dep_vertices = 0

    independent_vertices = []

    # initialize vertice containing map between vertices and indep index
    independent_vertex_indices = list(range(0, len(vertices)))

    for index, vertex in enumerate(vertices):
        independent_vertex = independent_vertex_position(vertex)

        if equal_vertices(vertex, independent_vertex):
            independent_vertices.append(vertex)
            continue
        else:

            ind_index = find_vertex_in_list(independent_vertices, independent_vertex)
            independent_vertex_indices[index] = ind_index

            num_indep_vertices -= 1
            num_dep_vertices += 1

    return independent_vertices, independent_vertex_indices


def create_custom_meshing_file(resolution):
    custom_name = 'custom_meshing_file.json'
    original_name = 'refined-meshing_opts.json'
    coarsening = 2

#    print(script_directory)
    original_name = script_directory + "/" + original_name

    with open(original_name) as original_file:
        meshing_opts = json.load(original_file)

        meshing_opts['marchingSquaresGridSize'] = int(2 ** coarsening * resolution)
        meshing_opts['marchingSquaresCoarsening'] = int(coarsening)

        with open(custom_name, 'w') as outfile:
            json.dump(meshing_opts, outfile)

    return custom_name


def generate_position_parameters(vertices):
    parameters = []

    for vertex in vertices:

        # first, normalize vertex
        normalized_vertex = [(vertex[0] + 1.0) / 2, (vertex[1] + 1.0) / 2]

        # then, verify if each coordinate is free and, if so, compute parameter
        if (not is_zero(normalized_vertex[0])) and (not is_zero(normalized_vertex[0] - 1.0)):
            parameters.append(normalized_vertex[0])

        if (not is_zero(normalized_vertex[1])) and (not is_zero(normalized_vertex[1] - 1.0)):
            parameters.append(normalized_vertex[1])

        # cause z is always at 0 coordinate
        #parameters.append(0.5)

    return parameters


def generate_parameters(independent_vertices, default_thickness, default_blending, custom_thickness_pairs):
    parameters = generate_position_parameters(independent_vertices)

    for i in range(0, len(independent_vertices)):
        parameters.append(default_thickness)

    for i in range(0, len(independent_vertices)):
        parameters.append(default_blending)

    offset_thickeness = len(parameters) - 2 * len(independent_vertices)

    # Now, for the extra parameters, find vertices and apply customized thickness
    for thickness_pair in custom_thickness_pairs:
        custom_thickness = thickness_pair[1]
        customized_nodes = thickness_pair[0]

        for node in customized_nodes:
            independent_vertex = independent_vertex_position(node)
            index = find_vertex_in_list(independent_vertices, independent_vertex)
            if index >= 0:
                parameters[offset_thickeness + index] = custom_thickness
            else:
                print("[Warning] vertex with customized thickness not found")

    return parameters


def generate_vertices_parameters(vertices, vertices_thickness, vertices_bending, custom_thickness_pairs):
    independent_vertices, independent_vertex_indices = extract_independent_vertices(vertices)

    parameters = generate_parameters(independent_vertices, vertices_thickness, vertices_bending, custom_thickness_pairs)

    return parameters

def inflate_hexagonal_box_smarter(input_path, vertices_thickness, vertices_bending, out_path,
                                  custom_thickness_pairs=[], resolution=64):
    # discover vertices
    vertices = []
    floats_pattern = re.compile(r'\-*\d+\.\d+')  # Compile a pattern to capture float values
    with open(input_path, 'r') as input_wire:
        wire_content = input_wire.readlines()
        for line in wire_content:

            if line.startswith('v '):
                floats = []
                for element in line.split():
                    try:
                        floats.append(float(element))
                    except:
                        continue

                vertices.append([floats[0], floats[1]])

    independent_vertices, independent_vertex_indices = extract_independent_vertices(vertices)

    parameters = generate_parameters(independent_vertices, vertices_thickness, vertices_bending, custom_thickness_pairs)

    parameters_string = ' '.join(str(param) for param in parameters)
    print(parameters_string)

    parameters_file_path = os.path.splitext(input_path)[0] + '.param'
    parameters_file = open(parameters_file_path, "w")
    parameters_file.write(parameters_string)
    parameters_file.close()

    custom_meshing_path = create_custom_meshing_file(resolution)

    cwd = script_directory
    cmd = [cwd + '/../../isosurface_inflator/isosurface_cli', '2D_doubly_periodic', input_path, '--paramsFile',
           parameters_file_path, '-m', custom_meshing_path, '--cheapPostprocessing' ,
           '-D', 'inflated.msh', '-R', 'replicated.msh', out_path]
#    (print cmd)
    check_output(cmd)

    if os.path.isfile(out_path):
        open(out_path + '.done', 'a').close()


def extract_vertices_from_edges_descriptions(edges_descriptions):
    extracted_vertices = []

    for index, edge_description in enumerate(edges_descriptions):
        vertex1 = edge_description[0]
        vertex2 = edge_description[1]

        position1 = find_vertex_in_list(extracted_vertices, vertex1)
        if position1 == -1:
            extracted_vertices.append(vertex1)

        position2 = find_vertex_in_list(extracted_vertices, vertex2)
        if position2 == -1:
            extracted_vertices.append(vertex2)

    return extracted_vertices


def remove_vertices_from_list(removed_vertices, original_vertices):
    vertices = list(original_vertices)

    for vertex in removed_vertices:
        index = find_vertex_in_list(vertices, vertex)
        if index >= 0:
            vertices.pop(index)

    return vertices


def connect_all_vertices(vertices):
    edges = []

    for index in range(0, len(vertices) - 1):
        edges.append([index, index + 1])

    return edges


def polygon_centroid(polygon_vertices):
    centroid = np.array([0.0, 0.0])

    for vertex in polygon_vertices:
        centroid[0] += vertex[0]
        centroid[1] += vertex[1]

    centroid[0] /= len(polygon_vertices)
    centroid[1] /= len(polygon_vertices)

    return centroid


def min_distance_point_line(point, line):
    x0 = point[0]
    y0 = point[1]
    x1 = line[0][0]
    y1 = line[0][1]
    x2 = line[1][0]
    y2 = line[1][1]
    line_norm = math.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1))

    distance = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / line_norm

    return distance


def add_new_vertices(new_vertices, vertices):
    for vertex in new_vertices:
        position = find_vertex_in_list(vertices, vertex)
        if position == -1:
            vertices.append(vertex)


def simplex_vertices_to_whole_parallelogram(simplex_vertices, parallelogram_side):
    l = parallelogram_side / 2
    vertices = list(simplex_vertices)

    # Operation 1: reflect against 30 degrees line
    vertices = np.array(vertices)
    reflected_vertices = reflect(30, vertices)

    vertices = vertices.tolist()
    add_new_vertices(reflected_vertices, vertices)
    vertices_step_1 = list(vertices)

    # Operation 2: reflect against y axis at x = l
    vertices -= np.array([l, 0])
    reflected_vertices = reflect(90, vertices)

    vertices = vertices.tolist()
    add_new_vertices(reflected_vertices, vertices)
    vertices += np.array([l, 0])

    # Operation 3: rotate 240 degrees and plug on top of current structure
    translated_vertices = vertices_step_1 - np.array([l, math.sqrt(3) / 3])
    rotated_vertices = rotate(240, translated_vertices)

    vertices = vertices.tolist()
    rotated_vertices += np.array([l, math.sqrt(3) / 3])
    add_new_vertices(rotated_vertices, vertices)

    # Operation 4:
    translated_vertices = vertices - np.array([2 * l, 0])
    reflected_vertices = reflect(120, translated_vertices)

    reflected_vertices += np.array([2 * l, 0])
    add_new_vertices(reflected_vertices, vertices)

    # Operation 5: transform parallelogram to square
    vertices = np.array(vertices)
    transformation_matrix = np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
    resulting_vertices = transformation_matrix * vertices.transpose()
    vertices = np.asarray(resulting_vertices.transpose())

    # Finally, transpose to origin
    vertices -= l

    return vertices


def simplex_to_whole_parallelogram(simplex_vertices, simplex_edges, parallelogram_side):
    l = parallelogram_side / 2
    vertices = list(simplex_vertices)
    edges = list(simplex_edges)

    # Operation 1: reflect against 30 degrees line
    vertices = np.array(vertices)
    reflected_vertices = reflect(30, vertices)

    reflected_edges = list(edges)
    vertices = vertices.tolist()
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges)
    vertices_step_1 = list(vertices)
    edges_step_1 = list(edges)

    # Operation 2: reflect against y axis at x = l
    vertices -= np.array([l, 0])
    reflected_vertices = reflect(90, vertices)

    reflected_edges = list(edges)
    vertices = vertices.tolist()
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges)
    vertices += np.array([l, 0])

    # Operation 3: rotate 240 degrees and plug on top of current structure
    translated_vertices = vertices_step_1 - np.array([l, math.sqrt(3) / 3])
    rotated_vertices = rotate(240, translated_vertices)

    rotated_edges = list(edges_step_1)
    vertices = vertices.tolist()
    rotated_vertices += np.array([l, math.sqrt(3) / 3])
    add_new_vertices_and_edges(rotated_vertices, rotated_edges, vertices, edges)

    # Operation 4:
    translated_vertices = vertices - np.array([2 * l, 0])
    reflected_vertices = reflect(120, translated_vertices)

    reflected_edges = list(edges)
    reflected_vertices += np.array([2 * l, 0])
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges)

    # Operation 5: transform parallelogram to square
    vertices = np.array(vertices)
    transformation_matrix = np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
    resulting_vertices = transformation_matrix * vertices.transpose()
    vertices = np.asarray(resulting_vertices.transpose())

    # Finally, transpose to origin
    vertices -= l

    return vertices, edges


def simplex_vertices_to_whole_square(simplex_vertices):
    vertices = list(simplex_vertices)

    # Operation 1: reflect against 45 degrees line
    vertices = np.array(vertices)
    reflected_vertices = reflect(45, vertices)

    vertices = vertices.tolist()
    add_new_vertices(reflected_vertices, vertices)

    # Operation 2: reflect against y axis at x = 0
    reflected_vertices = reflect(90, vertices)

    add_new_vertices(reflected_vertices, vertices)

    # Operation 3: reflect against x axis at y = 0
    reflected_vertices = reflect(0, vertices)

    add_new_vertices(reflected_vertices, vertices)

    return vertices


def simplex_polygon_to_whole_parallelogram(simplex_polygon, parallelogram_side):
    l = parallelogram_side / 2
    vertices = list(simplex_polygon)
    polygons = [simplex_polygon]

    # Operation 1: reflect against 30 degrees line
    vertices = np.array(vertices)
    reflected_vertices = reflect(30, vertices)
    polygons.append(reflected_vertices)  # polygons now has 2 pols
    polygons_step_1 = list(polygons)

    # Operation 2: reflect against y axis at x = l
    polygons_copy = list(polygons)
    for pol in polygons_copy:
        pol_vertices = pol - np.array([l, 0])
        reflected_vertices = reflect(90, pol_vertices)
        reflected_vertices += np.array([l, 0])
        polygons.append(reflected_vertices)
    # polygons now has 4 pols

    # Operation 3: rotate 240 degrees and plug on top of current structure
    polygons_copy = list(polygons_step_1)
    for pol in polygons_copy:
        translated_vertices = pol - np.array([l, math.sqrt(3) / 3])
        rotated_vertices = rotate(240, translated_vertices)
        rotated_vertices += np.array([l, math.sqrt(3) / 3])
        polygons.append(rotated_vertices)
    # polygons now has 6 pols

    # Operation 4:
    polygons_copy = list(polygons)
    for pol in polygons_copy:
        translated_vertices = pol - np.array([2 * l, 0])
        reflected_vertices = reflect(120, translated_vertices)
        reflected_vertices += np.array([2 * l, 0])
        polygons.append(reflected_vertices)
    # polygons now has 12 pols

    # Operation 5: transform parallelogram to square
    transformation_matrix = np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
    for index, pol in enumerate(polygons):
        vertices = np.array(pol)
        resulting_vertices = transformation_matrix * vertices.transpose()
        polygons[index] = np.asarray(resulting_vertices.transpose() - l)

    return polygons


def simplex_to_whole_square(simplex_vertices, simplex_edges):
    vertices = list(simplex_vertices)
    edges = list(simplex_edges)

    # Operation 1: reflect against 45 degrees line
    vertices = np.array(vertices)
    reflected_vertices = reflect(45, vertices)

    reflected_edges = list(edges)
    vertices = vertices.tolist()
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges)

    # Operation 2: reflect against y axis at x = 0
    reflected_vertices = reflect(90, vertices)

    reflected_edges = list(edges)
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges)

    # Operation 3: reflect against x axis at y = 0
    reflected_vertices = reflect(0, vertices)

    reflected_edges = list(edges)
    add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges)

    return vertices, edges


def rotate(theta_degrees, vertices):
    theta = math.radians(theta_degrees)
    R = np.matrix([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])

    reflected_vertices = R * np.array(vertices).transpose()

    return np.asarray(reflected_vertices.transpose())


def reflect(theta_degrees, vertices):
    theta = math.radians(theta_degrees)
    Rl = np.matrix([[math.cos(2 * theta), math.sin(2 * theta)], [math.sin(2 * theta), -math.cos(2 * theta)]])

    np_vertices = np.array(vertices)
    reflected_vertices = Rl * np_vertices.transpose()

    return np.asarray(reflected_vertices.transpose())


def add_new_vertices_and_edges(reflected_vertices, reflected_edges, vertices, edges):
    new_edges_description = []
    for edge in reflected_edges:
        v1 = edge[0]
        v2 = edge[1]

        reflected_vertex1 = reflected_vertices[v1]
        reflected_vertex2 = reflected_vertices[v2]

        new_edges_description.append([reflected_vertex1, reflected_vertex2])

    add_new_edges(new_edges_description, vertices, edges)


def triangle_incenter(triangle):
    p_a = np.array(triangle[0])
    p_b = np.array(triangle[1])
    p_c = np.array(triangle[2])

    a = np.linalg.norm(p_b - p_c)
    b = np.linalg.norm(p_a - p_c)
    c = np.linalg.norm(p_a - p_b)

    incenter = (a * p_a + b * p_b + c * p_c) / (a + b + c)

    return incenter


def polygon_to_edges_descriptions(polygon):
    edges_descriptions = []

    for i in range(0, len(polygon)):
        p1 = polygon[i]
        p2 = polygon[(i + 1) % len(polygon)]

        edges_descriptions.append([p1, p2])
    return edges_descriptions


def theoretical_rectangle(E_base, nu_base, volume_fraction):
    vertices = []
    k_base = E_base / (2 * (1 - nu_base))
    mu_base = E_base / (2 * (1 + nu_base))

    rho = volume_fraction
    k_lower = 0.0
    mu_lower = 0.0

    k_upper = k_base + (1 - rho) / (rho / (k_base + mu_base) - 1 / k_base)
    mu_upper = mu_base + (1 - rho) / (rho * (k_base + 2 * mu_base) / (2 * mu_base * (k_base + mu_base)) - 1 / mu_base)

    vertices.append([k_lower, mu_lower])
    vertices.append([k_lower, mu_upper])
    vertices.append([k_upper, mu_upper])
    vertices.append([k_upper, mu_lower])

    return vertices


def theoretical_triangle(E_base, nu_base, volume_fraction):
    rectangle_vertices = theoretical_rectangle(E_base, nu_base, volume_fraction)

    upper_vertex = rectangle_vertices[2]
    k = upper_vertex[0]
    mu = upper_vertex[1]

    # top vertex
    nu_top = (k - mu) / (k + mu)
    E_top = 2 * mu * (1 + nu_top)
    top_vertex = [nu_top, E_top]

    # left vertex
    left_vertex = [-1.0, 0.0]

    # right vertex
    right_vertex = [1.0, 0.0]

    return [left_vertex, top_vertex, right_vertex]


def add_polygons_incenters(polygons, vertices, edges):
    incenters_thickness_pairs = []

    for index, pol in enumerate(polygons):
        incenter = triangle_incenter(pol)
        thickness = min(min_distance_point_line(incenter, [pol[0], pol[1]]),
                        min_distance_point_line(incenter, [pol[1], pol[2]]),
                        min_distance_point_line(incenter, [pol[2], pol[0]]))

        incenters_thickness_pairs.append([[incenter], thickness])

        edges_descriptions = []
        for i in range(0, len(pol)):
            a = list(pol[i])
            b = incenter
            edges_descriptions.append([a, b])

        add_new_edges(edges_descriptions, vertices, edges)

    return incenters_thickness_pairs


def det_2D(a, b):
    return a[0] * b[1] - a[1] * b[0]


def edge_intersection(edge_1, edge_2):
    x_delta = (edge_1[0][0] - edge_1[1][0], edge_2[0][0] - edge_2[1][0])
    y_delta = (edge_1[0][1] - edge_1[1][1], edge_2[0][1] - edge_2[1][1])

    div = det_2D(x_delta, y_delta)
    if div == 0:
        raise Exception('lines do not intersect')

    d = (det_2D(*edge_1), det_2D(*edge_2))
    x = det_2D(d, x_delta) / div
    y = det_2D(d, y_delta) / div
    return [x, y]


def compute_volume_info(experiment_type, p1, p2, p3, p4):
    if experiment_type == "negative":
        height = 1 * math.sqrt(3) / 3.0
        s = 2 * height

        volume_triangle = math.sqrt(3) / 4 * (p1 * s) ** 2
        volume_pillars = p1 * p3 * p4 * s * (1 - math.sqrt(3) / 2 * s * p1)
        volume = volume_triangle + volume_pillars

        total_volume = math.sqrt(3) / 4 * s ** 2

        volume_fraction = volume / total_volume

    else:
        r = p1 * math.sqrt(3) / 3
        x = p1

        volume_triangle = r * p1
        volume_pillars = 2 * math.sqrt(3) / 3 * (1 - p1) * p1 * p3 * p4

        volume_fraction = p1 ** 2 + 2 * (1 - p1) * p1 * p3 * p4

        volume = volume_triangle + volume_pillars
        total_volume = math.sqrt(3) / 3
        assert abs(volume_fraction - volume / total_volume) < Tolerance

    print("Volume is: ", volume)
    print("Volume fraction is ", volume_fraction)

    return volume, volume_fraction


def generate_positive_poisson_topology_and_thickness_info(p1, p2, p4):
    triangle_side_ratio = p1
    num_pillars = p2
    pillar_area_ratio = 1.0
    thickness_ratio = p4
    parallelogram_side = 2.0
    l = parallelogram_side / 2.0

    triangle_side = triangle_side_ratio * 2
    thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars) / 2  # use the radius
    # graph constructed on left triangle
    vertices = []
    edges = []

    # define important vertices of parallelogram
    origin = np.array([0, 0])
    triangle_centroid = np.array([l, l * math.sqrt(3) / 3.0])

    # First side: create pillars of the left side (counter-clock orientation)
    #     /.
    #    /  .
    #   /.....
    parallelogram_start = np.array([(l + triangle_side / 2.0) / 2.0, (l + triangle_side / 2.0) * math.sqrt(3) / 2.0])
    parallelogram_end = np.array([(l - triangle_side / 2.0) / 2.0, (l - triangle_side / 2.0) * math.sqrt(3) / 2.0])

    radius = triangle_side * math.sqrt(3) / 6
    f = triangle_centroid - [0, radius]
    triangle_start = f + [0, math.sqrt(3) * triangle_side / 2.0]
    triangle_end = f - [triangle_side / 2.0, 0]

    unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
    unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

    parallelogram_start_offset = parallelogram_start + thickness * unit_vector
    parallelogram_end_offset = parallelogram_end - thickness * unit_vector
    triangle_start_offset = triangle_start + thickness * unit_vector
    triangle_end_offset = triangle_end - thickness * unit_vector

    # add new vertices and edges to current sets
    new_edges_description = create_pillars([parallelogram_start_offset, parallelogram_end_offset],
                                                  [triangle_start_offset, triangle_end_offset], num_pillars)

    # save pillar nodes
    pillar_nodes = []
    for edge_description in new_edges_description:
        pillar_nodes.append(edge_description[0])
        pillar_nodes.append(edge_description[1])

    # new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
    add_new_edges(new_edges_description, vertices, edges)

    # Second side: create pillars of the bottom side (counter-clock orientation)
    #     ..
    #    .  .
    #   ______
    parallelogram_start = np.array([l - triangle_side / 2.0, 0])
    parallelogram_end = np.array([l + triangle_side / 2.0, 0])

    triangle_start = f - [triangle_side / 2.0, 0]
    triangle_end = f + [triangle_side / 2.0, 0]

    pillar_example = [[parallelogram_start, triangle_start]]

    unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
    unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

    parallelogram_start_offset = parallelogram_start + thickness * unit_vector
    parallelogram_end_offset = parallelogram_end - thickness * unit_vector
    triangle_start_offset = triangle_start + thickness * unit_vector
    triangle_end_offset = triangle_end - thickness * unit_vector

    # add new vertices and edges to current sets
    new_edges_description = create_pillars([parallelogram_start_offset, parallelogram_end_offset],
                                                  [triangle_start_offset, triangle_end_offset],
                                                  num_pillars)

    # save pillar nodes
    for edge_description in new_edges_description:
        pillar_nodes.append(edge_description[0])
        pillar_nodes.append(edge_description[1])

    # new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
    add_new_edges(new_edges_description, vertices, edges)

    # Third side: create pillars of the right side (counter-clock orientation)
    #     .\
    #    .  \
    #   .....\
    parallelogram_start = np.array(
        [2 * l - (l - triangle_side / 2.0) / 2, (l - triangle_side / 2.0) * math.sqrt(3) / 2])
    parallelogram_end = np.array([2 * l - (l + triangle_side / 2.0) / 2, (l + triangle_side / 2.0) * math.sqrt(3) / 2])

    triangle_start = f + [triangle_side / 2.0, 0]
    triangle_end = f + [0, math.sqrt(3) * triangle_side / 2.0]

    unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
    unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

    parallelogram_start_offset = parallelogram_start + thickness * unit_vector
    parallelogram_end_offset = parallelogram_end - thickness * unit_vector
    triangle_start_offset = triangle_start + thickness * unit_vector
    triangle_end_offset = triangle_end - thickness * unit_vector

    # add new vertices and edges to current sets
    new_edges_description = create_pillars([parallelogram_start_offset, parallelogram_end_offset],
                                                  [triangle_start_offset, triangle_end_offset],
                                                  num_pillars)

    # save nodes on hypotenuse
    hypotenuse_nodes = []
    for edge_description in new_edges_description:
        hypotenuse_nodes.append(edge_description[0])
        hypotenuse_nodes.append(edge_description[1])

    # new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
    add_new_edges(new_edges_description, vertices, edges)

    # Now, adding edges on the triangle
    triangle = [f + [0, math.sqrt(3) * triangle_side / 2.0], f - [triangle_side / 2.0, 0], f + [triangle_side / 2.0, 0]]
    new_edges_description = create_triangle_edges(triangle, num_pillars, thickness, thickness)
    add_new_edges(new_edges_description, vertices, edges)

    # Now, transform to square every vertex we have:
    vertices = np.array(vertices)
    transformation_matrix = np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
    resulting_vertices = transformation_matrix * vertices.transpose()
    vertices = np.asarray(resulting_vertices.transpose())

    # Finally, transpose to origin and reflect through the diagonal 'y = -x'
    vertices -= 1.0

    reflected_vertices = vertices[:, [1, 0]] * (-1)

    # saving vertices on hypotenuses
    hypotenuse_nodes = np.array(hypotenuse_nodes)
    resulting_hypotenuse_nodes = transformation_matrix * hypotenuse_nodes.transpose()
    hypotenuse_nodes = np.asarray(resulting_hypotenuse_nodes.transpose())
    hypotenuse_nodes -= 1.0
    reflected_hypotenuse_nodes = hypotenuse_nodes[:, [1, 0]] * (-1)
    hypotenuse_nodes = hypotenuse_nodes.tolist() + reflected_hypotenuse_nodes.tolist()

    # saving vertices on pillars
    pillar_nodes = np.array(pillar_nodes)
    resulting_pillar_nodes = transformation_matrix * pillar_nodes.transpose()
    pillar_nodes = np.asarray(resulting_pillar_nodes.transpose())
    pillar_nodes -= 1.0
    reflected_pillar_nodes = pillar_nodes[:, [1, 0]] * (-1)
    pillar_nodes = pillar_nodes.tolist() + reflected_pillar_nodes.tolist()

    # transforming pillar example
    pillar_example = np.array(pillar_example)
    resulting_pillar_example = transformation_matrix * pillar_example.transpose()
    pillar_example = np.asarray(resulting_pillar_example.transpose())
    delta_y = abs(pillar_example[1][1] - pillar_example[0][1])
    norm_example = np.linalg.norm(pillar_example[1] - pillar_example[0])
    thickness_correction_factor = delta_y / norm_example
    # print("Thickness correction factor: " + str(thickness_correction_factor))

    new_edges_description = []
    for edge in edges:
        v1 = edge[0]
        v2 = edge[1]

        reflected_vertex1 = reflected_vertices[v1]
        reflected_vertex2 = reflected_vertices[v2]

        new_edges_description.append([reflected_vertex1, reflected_vertex2])

    vertices = vertices.tolist()
    add_new_edges(new_edges_description, vertices, edges)

    # triangle is also transformed and reflected
    triangle = np.array(triangle)
    resulting_triangle = transformation_matrix * triangle.transpose()
    transformed_triangle = np.asarray(resulting_triangle.transpose())

    transformed_triangle -= 1.0
    reflected_triangle = transformed_triangle[:, [1, 0]] * (-1)

    incenter_triangle_pairs = add_polygons_incenters([transformed_triangle, reflected_triangle], vertices, edges)

    custom_thickness_pairs = incenter_triangle_pairs + [[hypotenuse_nodes, float(thickness) * math.sqrt(2)],
                               [pillar_nodes, float(thickness) * thickness_correction_factor]]

    return vertices, edges, custom_thickness_pairs



