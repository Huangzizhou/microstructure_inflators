import math
import os
import re
from subprocess import call

import numpy as np

Tolerance = 1e-3


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


def create_triangle_edges(triangle, n):
    edges_descriptions = []

    for i in range(0, 3):
        segment_start = triangle[i % 3]
        segment_end = triangle[(i + 1) % 3]
        create_subdivided_segment(segment_start, segment_end, n, edges_descriptions)

    return edges_descriptions


def create_subdivided_segment(segment_start, segment_end, n, edges_descriptions):
    vector = segment_end - segment_start
    last_point = segment_start
    for t in range(1, n):
        new_point = segment_start + 1.0 * t / (n - 1) * vector

        edges_descriptions.append([last_point, new_point])
        last_point = new_point


def create_subdivided_segment_with_constant_spacing(segment_start, segment_end, n, spacing, edges_descriptions):
    vector = (segment_end - segment_start)
    direction = vector / np.linalg.norm(vector)
    last_point = segment_start
    for t in range(1, n-1):
        new_point = segment_start + spacing * t * direction

        edges_descriptions.append([last_point, new_point])
        last_point = new_point

    edges_descriptions.append([last_point, segment_end])


def create_pillars(segment1, segment2, n):
    edges_descriptions = []

    # vectors representing each segment
    vector1 = segment1[1] - segment1[0]
    vector2 = segment2[1] - segment2[0]

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


def find_vertex_in_list(vertices, v):
    new_x = v[0]
    new_y = v[1]

    vertex_position = -1

    for index, vertex in enumerate(vertices):
        x = vertex[0]
        y = vertex[1]

        if (abs(x - new_x) < Tolerance) and (abs(y - new_y) < Tolerance):
            vertex_position = index

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
        cmd = [cwd + '/../../isosurface_inflator/default_parameters', '2D_triply_periodic', input_path]
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
            print "Warning: no vertex in this position", node
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

    cmd = [cwd + '/../../isosurface_inflator/isosurface_cli', '2D_triply_periodic', input_path, '--params',
           parameters_string, '-m', 'refined-meshing_opts.json', '-D', 'inflated.msh', '-R', 'replicated.msh', out_path]
    # print cmd
    call(cmd)


def inflate_hexagonal_box_smarter(input_path, vertices_thickness, vertices_bending, out_path, custom_thickness_pairs = []):
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

    # run default_parameters binary to obtain default parameters
    cwd = os.getcwd()
    with open('default-parameters.txt', 'w') as out_log:
        cmd = [cwd + '/../../isosurface_inflator/default_parameters', '2D_triply_periodic', input_path]
        call(cmd, stdout=out_log)

    with open('default-parameters.txt', 'r') as out_log:
        # read outlog to obtain default parameters and modify them
        content = out_log.readlines()
        parameters = []
        for line_index in range(1, len(content)):
            for element in content[line_index].split():
                # print element
                parameters.append(float(element))

    offset_thickeness = len(parameters) - 2 * len(vertices)
    for i in range(offset_thickeness, offset_thickeness + len(vertices)):
        parameters[i] = vertices_thickness

    offset_bending = len(parameters) - len(vertices)
    for i in range(offset_bending, offset_bending + len(vertices)):
        parameters[i] = vertices_bending

    # Now, for the extra parameters, find vertices and apply customized thickness
    for thickness_pair in custom_thickness_pairs:
        custom_thickness = thickness_pair[1]
        customized_nodes = thickness_pair[0]

        for node in customized_nodes:
            index = find_vertex_in_list(vertices, node)
            if index >= 0:
                parameters[offset_thickeness + index] = custom_thickness
            else:
                print "[Warning] vertex with customized thickness not found"

    parameters_string = ' '.join(str(param) for param in parameters)

    parameters_file_path = os.path.splitext(input_path)[0] + '.param'
    parameters_file = open(parameters_file_path, "w")
    parameters_file.write(parameters_string)
    parameters_file.close()

    cmd = [cwd + '/../../isosurface_inflator/isosurface_cli', '2D_triply_periodic', input_path, '--params',
           parameters_string, '-m', 'refined-meshing_opts.json', '-D', 'inflated.msh', '-R', 'replicated.msh', out_path]
    # print cmd
    call(cmd)


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
        edges.append([index, index+1])

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
    line_norm = math.sqrt((y2-y1)*(y2-y1) + (x2-x1)*(x2-x1))

    distance = abs((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1) / line_norm

    return distance


def add_new_vertices(new_vertices, vertices):

    for vertex in new_vertices:
        position = find_vertex_in_list(vertices, vertex)
        if position == -1:
            vertices.append(vertex)


def simplex_vertices_to_whole_parallelogram(simplex_vertices, parallelogram_side):
    l = parallelogram_side/2
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
    polygons.append(reflected_vertices) # polygons now has 2 pols
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

    reflected_vertices = R * vertices.transpose()

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

    incenter = (a*p_a + b*p_b + c*p_c) / (a + b + c)

    return incenter


def polygon_to_edges_descriptions(polygon):
    edges_descriptions = []

    for i in range(0, len(polygon)):
        p1 = polygon[i]
        p2 = polygon[(i+1) % len(polygon)]

        edges_descriptions.append([p1, p2])
    return edges_descriptions