#!/usr/bin/env python
import sys
import math

import numpy as np

import hexlib


def get_thickness(min_thickness_ratio, max_thickness_ratio, index, num_pillars):
    w_min = min_thickness_ratio * (pillar_area / num_pillars)
    w_max = max_thickness_ratio * (pillar_area / num_pillars)

    if num_pillars < 3:
        return w_min
    elif num_pillars % 2 == 1:
        # odd
        m = (num_pillars + 1) / 2

        angular_coeff = (w_max - w_min) / (m - 1)
        linear_coeff = w_min
    else:
        # even
        m = num_pillars / 2
        angular_coeff = (w_max - w_min) / (m - 1)
        linear_coeff = w_min

    if index > m:
        result = linear_coeff + (num_pillars - index) * angular_coeff
    else:
        result = linear_coeff + (index - 1) * angular_coeff

    return result

def get_spacing(min_thickness, max_thickness, num_pillars, total_pillar_region):
    w_min = min_thickness_ratio * (total_pillar_region / num_pillars)
    w_max = max_thickness_ratio * (total_pillar_region / num_pillars)

    if num_pillars < 3:
        spacing = (total_pillar_region - w_min) / (num_pillars - 1) - w_min
        print "Thickness: " + str(w_min)
        print "Spacing: " + str(spacing)
        return spacing
    if num_pillars % 2 == 1:
        # odd
        m = (num_pillars + 1) / 2
        angular_coeff = (w_max - w_min) / (m - 1)
        linear_coeff = w_min

        sum = (2*m - 1) * linear_coeff + angular_coeff * (1 -2*m + m ** 2)
    else:
        # even
        m = num_pillars / 2
        angular_coeff = (w_max - w_min) / (m - 1)
        linear_coeff = w_min

        sum = 2*m * linear_coeff + angular_coeff * (-m + m ** 2)

    # computes constant spacing
    result = (total_pillar_region - sum) / (num_pillars - 1)

    return result

def translate_at(origin, offset, vertices, edges, edges_descriptions):
    translated_vertices = vertices - origin + offset
    for edge in edges:
        v1 = translated_vertices[edge[0]]
        v2 = translated_vertices[edge[1]]

        edges_descriptions.append([v1, v2])

    return translated_vertices


def rotate_at(center_of_rotation, degrees, vertices, edges, edges_descriptions):
    translated_vertices = vertices - center_of_rotation
    rotated_vertices = hexlib.rotate(degrees, translated_vertices) + center_of_rotation
    for edge in edges:
        v1 = rotated_vertices[edge[0]]
        v2 = rotated_vertices[edge[1]]

        edges_descriptions.append([v1, v2])

    return rotated_vertices


def create_pillar_edge(edge):
    crossing_edge = [[-1.0, 0.0], [1.0, 0.0]]
    edge_intersection_point = hexlib.edge_intersection(edge, crossing_edge)
    resulting_edge = [edge[0], edge_intersection_point]
    return resulting_edge


def create_pillars_with_constant_spacing_and_thickness(line1, line2, num_pillars, min_thickness_ratio, max_thickness_ratio, pillar_area):
    edges_descriptions = []
    unit_vector = (line1[1] - line1[0]) / np.linalg.norm(line1[1] - line1[0])
    pillar_polygons = []

    # initial points for each side
    next_point1 = line1[0]
    next_point2 = line2[0]

    for t in range(0, num_pillars):
        thickness = get_thickness(min_thickness_ratio, max_thickness_ratio, t+1, num_pillars)
        gap = get_spacing(min_thickness_ratio, max_thickness_ratio, num_pillars, pillar_area)

        beginning_pillar1 = next_point1
        midpoint_pillar1 = beginning_pillar1 + thickness * unit_vector / 2
        end_pillar1 = beginning_pillar1 + thickness * unit_vector

        beginning_pillar2 = next_point2
        midpoint_pillar2 = beginning_pillar2 + thickness * unit_vector / 2
        end_pillar2 = beginning_pillar2 + thickness * unit_vector

        top_edge1 = [beginning_pillar1, midpoint_pillar1]
        top_edge2 = [midpoint_pillar1, end_pillar1]

        pillar_edge = [midpoint_pillar1, midpoint_pillar2]
        resulting_pillar_edge = create_pillar_edge(pillar_edge)

        pillar_edges = [top_edge1, top_edge2, resulting_pillar_edge]

        pillar_polygons.append([midpoint_pillar1, beginning_pillar1, beginning_pillar2, resulting_pillar_edge[1]])

        next_point1 = end_pillar1 + gap * unit_vector
        next_point2 = end_pillar2 + gap * unit_vector

        if t < num_pillars - 1:
            gap_edge = [end_pillar1, next_point1]
            pillar_edges.append(gap_edge)

        for pillar_edge in pillar_edges:
            if not pillar_edge:
                pass
            else:
                edges_descriptions.append(pillar_edge)

    return edges_descriptions, pillar_polygons



if len(sys.argv) != 9:
    print "usage: ./auxetic-gradient-creator.py <size of triangles> <num of pillars> <relative size of pillar area (chirality)> <min relative thickness (of pillars)> <max relative thickness (of pillars)> <ninja coefficient>" \
          "<output wire> <output mesh>"
    print "example: ./auxetic-gradient-creator.py 0.8 3 0.6 0.9 0.8 0.9 output.wire output.msh"
    sys.exit(-1)

#       .b
#      /\
#     /  \
#    /    \
# a /______\ c
#
triangle_side_ratio = float(sys.argv[1])
num_pillars = int(sys.argv[2])
pillar_area_ratio = float(sys.argv[3])
min_thickness_ratio = float(sys.argv[4])
max_thickness_ratio = float(sys.argv[5])
ninja_factor = float(sys.argv[6])
out_wire = sys.argv[7]
out_mesh = sys.argv[8]
parallelogram_side = 3.0
s = parallelogram_side / 3.0
p1 = float(sys.argv[1])
p2 = int(sys.argv[2])
p3 = float(sys.argv[3])
p4 = float(sys.argv[4])
p5 = float(sys.argv[5])
p6 = float(sys.argv[6])


print "Constructing " + out_wire + " ..."

# graph to be constructed
vertices = []
edges = []
triangles = []
pillar_polygons = []

# define important vertices of simplex used to build entire parallelogram structure
origin = np.array([0, 0])
a = origin.copy()
b = np.array([s/2.0, s * math.sqrt(3)/2.0])
c = np.array([s, 0])

# define vertices of triangle
triangle_side = triangle_side_ratio * s
triangle_y_position = s * math.sqrt(3)/2.0 * (1-triangle_side_ratio)

q1 = [s/2.0 * (1-triangle_side_ratio), triangle_y_position]
q2 = [s/2.0 * (1+triangle_side_ratio), triangle_y_position]
w  = [q2[0] - p1*p3*s , triangle_y_position]
ba_unit = (b - a) / np.linalg.norm((b-a))
z = ba_unit * p1*p3*s*(1-p6) + w

# reflect q1 and q2 against the x axis
q1_reflected = [q2[0], -q2[1]]
q2_reflected = [q1[0], -q1[1]]
w_reflected  = [q2_reflected[0] + p1*p3*s , -triangle_y_position]
ab_unit = (a-b) / np.linalg.norm(a-b)
z_reflected = ab_unit * p1*p3*s*(1-p6) + w_reflected

# computing thickness and spacing
pillar_area = np.linalg.norm(z - q2)

# CREATE VERTICES IN THE INITIAL (ORIGINAL) EQUILATERAL TRIANGLE
edges_descriptions, original_pillar_polygons = create_pillars_with_constant_spacing_and_thickness(np.array([q2, z]), np.array([z_reflected, q2_reflected]), num_pillars, min_thickness_ratio, max_thickness_ratio, pillar_area)
edges_descriptions.append([z, q1])
pillar_polygons += original_pillar_polygons
hexlib.add_new_edges(edges_descriptions, vertices, edges)

triangles.append([q1, b, z])
triangles.append([q2, b, z])

# CREATE BOTTOM-LEFT PART OF THE STRUCTURE
edges_descriptions = []
rotate_at(b, 60, vertices, edges, edges_descriptions)
rotate_at(b, 120, vertices, edges, edges_descriptions)

# do the same with hexagon vertices
hexagon_vertices = [q1, z, q2]
hexagon_edges = [[0, 1], [1, 2]]
hexagon_edges_descriptions = [[q1, z], [z, q2]]
rotate_at(b, 60, hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(b, 120, hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_vertices = hexlib.extract_vertices_from_edges_descriptions(hexagon_edges_descriptions)
for i in range(0, 6):
    triangles.append([b, hexagon_vertices[i], hexagon_vertices[i+1]])

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    rotated_60 = rotate_at(b, 60, pillar_polygon, [], [])
    rotated_120 = rotate_at(b, 120, pillar_polygon, [], [])
    pillar_polygons += [rotated_60, rotated_120]


# CREATE TOP-LEFT PART OF THE STRUCTURE
center = np.array([2*s, s*math.sqrt(3)])
translated_vertices = translate_at(b, center, vertices, edges, edges_descriptions)
rotate_at(center, 60, translated_vertices, edges, edges_descriptions)
rotate_at(center, 120, translated_vertices, edges, edges_descriptions)
rotate_at(center, 180, translated_vertices, edges, edges_descriptions)
rotate_at(center, 240, translated_vertices, edges, edges_descriptions)
rotate_at(center, 300, translated_vertices, edges, edges_descriptions)

# do the same with hexagon vertices
hexagon_vertices = [q1, z, q2]
hexagon_edges = [[0, 1], [1, 2]]
hexagon_edges_descriptions = []
translated_hexagon_vertices = translate_at(b, center, hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 60, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 120, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 180, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 240, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 300, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 360, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_vertices = hexlib.extract_vertices_from_edges_descriptions(hexagon_edges_descriptions)
for i in range(0, 11):
    triangles.append([center, hexagon_vertices[i], hexagon_vertices[i+1]])
triangles.append([center, hexagon_vertices[0], hexagon_vertices[11]])

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    translated_pillar_polygon = translate_at(b, center, pillar_polygon, [], [])
    rotated_60 = rotate_at(center, 60, translated_pillar_polygon, [], [])
    rotated_120 = rotate_at(center, 120, translated_pillar_polygon, [], [])
    rotated_180 = rotate_at(center, 180, translated_pillar_polygon, [], [])
    rotated_240 = rotate_at(center, 240, translated_pillar_polygon, [], [])
    rotated_300 = rotate_at(center, 300, translated_pillar_polygon, [], [])
    pillar_polygons += [translated_pillar_polygon, rotated_60, rotated_120, rotated_180, rotated_240, rotated_300]


# CREATE TOP-RIGHT PART OF THE STRUCTURE
center = np.array([3*s + s/2, 3.0/2*s*math.sqrt(3)])
translated_vertices = translate_at(b, center, vertices, edges, edges_descriptions)
rotate_at(center, 60, translated_vertices, edges, edges_descriptions)
rotate_at(center, -60, translated_vertices, edges, edges_descriptions)

# do the same with hexagon vertices
hexagon_vertices = [q1, z, q2]
hexagon_edges = [[0, 1], [1, 2]]
hexagon_edges_descriptions = []
translated_hexagon_vertices = translate_at(b, center, hexagon_vertices, hexagon_edges, [])
rotate_at(center, -60, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_edges_descriptions += [[translated_hexagon_vertices[0], translated_hexagon_vertices[1]], [translated_hexagon_vertices[1], translated_hexagon_vertices[2]]]
rotate_at(center, 60, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_vertices = hexlib.extract_vertices_from_edges_descriptions(hexagon_edges_descriptions)
for i in range(0, 6):
    triangles.append([center, hexagon_vertices[i], hexagon_vertices[i+1]])

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    translated_pillar_polygon = translate_at(b, center, pillar_polygon, [], [])
    rotated_60 = rotate_at(center, 60, translated_pillar_polygon, [], [])
    rotated_300 = rotate_at(center, -60, translated_pillar_polygon, [], [])
    pillar_polygons += [translated_pillar_polygon, rotated_60, rotated_300]


# CREATE BOTTOM PART OF THE STRUCTURE
center = np.array([2*s, 0.0])
translated_vertices = translate_at(b, center, vertices, edges, [])
rotate_at(center, 120, translated_vertices, edges, edges_descriptions)
rotate_at(center, 180, translated_vertices, edges, edges_descriptions)
rotate_at(center, 240, translated_vertices, edges, edges_descriptions)

# do the same with hexagon vertices
hexagon_vertices = [q1, z, q2]
hexagon_edges = [[0, 1], [1, 2]]
hexagon_edges_descriptions = []
translated_hexagon_vertices = translate_at(b, center, hexagon_vertices, hexagon_edges, [])
rotate_at(center, 120, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 180, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 240, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_vertices = hexlib.extract_vertices_from_edges_descriptions(hexagon_edges_descriptions)
for i in range(0, 6):
    triangles.append([center, hexagon_vertices[i], hexagon_vertices[i+1]])

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    translated_pillar_polygon = translate_at(b, center, pillar_polygon, [], [])
    rotated_120 = rotate_at(center, 120, translated_pillar_polygon, [], [])
    rotated_180 = rotate_at(center, 180, translated_pillar_polygon, [], [])
    rotated_240 = rotate_at(center, 240, translated_pillar_polygon, [], [])
    pillar_polygons += [rotated_120, rotated_180, rotated_240]


# CREATE RIGHT PART OF THE STRUCTURE
center = np.array([3*s + s/2, s*math.sqrt(3)/2])
translated_vertices = translate_at(b, center, vertices, edges, [])
rotate_at(center, 180, translated_vertices, edges, edges_descriptions)
rotate_at(center, 240, translated_vertices, edges, edges_descriptions)
rotate_at(center, 300, translated_vertices, edges, edges_descriptions)

# do the same with hexagon vertices
hexagon_vertices = [q1, z, q2]
hexagon_edges = [[0, 1], [1, 2]]
hexagon_edges_descriptions = []
translated_hexagon_vertices = translate_at(b, center, hexagon_vertices, hexagon_edges, [])
rotate_at(center, 180, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 240, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 300, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_vertices = hexlib.extract_vertices_from_edges_descriptions(hexagon_edges_descriptions)
for i in range(0, 6):
    triangles.append([center, hexagon_vertices[i], hexagon_vertices[i+1]])

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    translated_pillar_polygon = translate_at(b, center, pillar_polygon, [], [])
    rotated_180 = rotate_at(center, 180, translated_pillar_polygon, [], [])
    rotated_240 = rotate_at(center, 240, translated_pillar_polygon, [], [])
    rotated_300 = rotate_at(center, 300, translated_pillar_polygon, [], [])
    pillar_polygons += [rotated_180, rotated_240, rotated_300]

# ADD ALL EDGES DESCRIPTIONS
hexlib.add_new_edges(edges_descriptions, vertices, edges)

# TRANSFORM TO FINAL SQUARE SHAPE in [-1, 1] x [-1, 1] cell
vertices = np.array(vertices)
transformation_matrix = 2./3 * np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
#transformation_matrix = np.matrix('1.0 0.0; 0.0 1.0')
resulting_vertices = transformation_matrix * vertices.transpose()
vertices = np.asarray(resulting_vertices.transpose()) - 1.0

# deal with triangles
for index, triangle in enumerate(triangles):
    resulting_triangle = transformation_matrix * np.array(triangle).transpose()
    triangle = np.asarray(resulting_triangle.transpose()) - 1.0
    triangles[index] = triangle

vertices = list(vertices)
incenters_thickness_pairs = hexlib.add_polygons_incenters(triangles, vertices, edges)

# deal with pillars
pillar_nodes_custom_pairs = []
for index, pillar_polygon in enumerate(pillar_polygons):
    resulting_pillar_polygon = transformation_matrix * np.array(pillar_polygon).transpose()
    resulting_pillar_polygon = np.asarray(resulting_pillar_polygon.transpose()) - 1.0

    nodes_thickness = hexlib.min_distance_point_line(resulting_pillar_polygon[0], [resulting_pillar_polygon[1], resulting_pillar_polygon[2]])
    pillar_nodes = [resulting_pillar_polygon[0], resulting_pillar_polygon[3]]

    pillar_nodes_custom_pairs.append([pillar_nodes, nodes_thickness])


# FINALLY, CREATE WIRE
hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

thickness =  get_thickness(min_thickness_ratio, max_thickness_ratio, 1, num_pillars)
spacing =  get_spacing(min_thickness_ratio, max_thickness_ratio, num_pillars, pillar_area)

#TODO: find, in microstructure, smaller and larger pillar tickness, so we can compute thickness void and real thickness with more precision
thickness_void = (triangle_side*pillar_area_ratio - num_pillars*thickness) / (num_pillars - 1)
min_resolution = 2 * max(3 / thickness_void, 3 / thickness)
chosen_resolution = math.pow(2, math.ceil(math.log(min_resolution) / math.log(2)))

if chosen_resolution > 2048:
        print "Resolution of " + str(chosen_resolution) + "is too big"
        print "Skipping experiment!"
        exit()

if chosen_resolution < 64:
    chosen_resolution = 64

print "Thickness void: " + str(thickness_void)
print "Minimum resolution: " + str(min_resolution)
print "Chosen resolution: " + str(chosen_resolution)

hexlib.inflate_hexagonal_box_smarter(out_wire, 0.00001, 0.0, out_mesh, incenters_thickness_pairs + pillar_nodes_custom_pairs, chosen_resolution)
