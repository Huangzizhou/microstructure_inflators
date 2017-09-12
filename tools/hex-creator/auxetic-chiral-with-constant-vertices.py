#!/usr/bin/env python
import sys
import math

import numpy as np

import hexlib


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


def create_pillars_with_constant_spacing_and_thickness(line1, line2, pillar_area_ratio, num_pillars, thickness, gap):
    edges_descriptions = []
    unit_vector = (line1[1] - line1[0]) / np.linalg.norm(line1[1] - line1[0])
    pillar_polygons = []

    # initial points for each side
    next_point1 = line1[0]
    next_point2 = line2[1] + (line2[0] - line2[1]) * pillar_area_ratio

    for t in range(0, num_pillars):
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

    # need to add final edge connecting next_point1 to end of line1 (only if they differ)
    final_edge = [end_pillar1, line1[1]]
    edges_descriptions.append(final_edge)

    return edges_descriptions, pillar_polygons



if len(sys.argv) != 7:
    print "usage: ./auxetic-chiral-with-constant-vertices.py <size of triangles> <num of pillars> <relative size of pillar area (chirality)> <relative thickness (of pillars)> " \
          "<output wire> <output mesh>"
    print "example: ./auxetic-chiral-with-constant-vertices.py 0.8 3 0.6 0.9 output.wire output.msh"
    sys.exit(-1)

# . b
#     . '  |  ' .
# a._______|_______. c
#
triangle_side_ratio = float(sys.argv[1])
num_pillars = int(sys.argv[2])
pillar_area_ratio = float(sys.argv[3])
thickness_ratio = float(sys.argv[4])
out_wire = sys.argv[5]
out_mesh = sys.argv[6]
parallelogram_side = 3.0
s = parallelogram_side / 3.0

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

# reflect q1 and q2 against the x axis
q1_reflected = [q1[0], -q1[1]]
q2_reflected = [q2[0], -q2[1]]

thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars)
spacing = (pillar_area_ratio * triangle_side - thickness) / (num_pillars - 1) - thickness

edges_descriptions, original_pillar_polygons = create_pillars_with_constant_spacing_and_thickness(np.array([q2, q1]), np.array([q2_reflected, q1_reflected]), pillar_area_ratio, num_pillars, thickness, spacing)
pillar_polygons += original_pillar_polygons
hexlib.add_new_edges(edges_descriptions, vertices, edges)

edges_descriptions = []
rotate_at(b, 60, vertices, edges, edges_descriptions)
rotate_at(b, 120, vertices, edges, edges_descriptions)

# do the same with hexagon vertices
hexagon_vertices = [q1, q2]
hexagon_edges = [[0, 1]]
hexagon_edges_descriptions = [[q1, q2]]
rotate_at(b, 60, hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(b, 120, hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_vertices = hexlib.extract_vertices_from_edges_descriptions(hexagon_edges_descriptions)
triangles.append([hexagon_vertices[0], hexagon_vertices[1], hexagon_vertices[2]])
triangles.append([hexagon_vertices[0], hexagon_vertices[2], hexagon_vertices[3]])
triangles_to_be_reflected = [0, 1]

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    rotated_60 = rotate_at(b, 60, pillar_polygon, [], [])
    rotated_120 = rotate_at(b, 120, pillar_polygon, [], [])
    pillar_polygons += [rotated_60, rotated_120]


center = np.array([2*s, s*math.sqrt(3)])
translated_vertices = translate_at(b, center, vertices, edges, edges_descriptions)
rotate_at(center, 60, translated_vertices, edges, edges_descriptions)
rotate_at(center, 120, translated_vertices, edges, edges_descriptions)
rotate_at(center, 180, translated_vertices, edges, edges_descriptions)
rotate_at(center, 240, translated_vertices, edges, edges_descriptions)
rotate_at(center, 300, translated_vertices, edges, edges_descriptions)

# do the same with hexagon vertices
hexagon_vertices = [q1, q2]
hexagon_edges = [[0, 1]]
hexagon_edges_descriptions = []
translated_hexagon_vertices = translate_at(b, center, hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 60, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 120, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 180, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 240, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 300, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 360, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_vertices = hexlib.extract_vertices_from_edges_descriptions(hexagon_edges_descriptions)
triangles.append([hexagon_vertices[0], hexagon_vertices[1], hexagon_vertices[2]])
triangles.append([hexagon_vertices[3], hexagon_vertices[4], hexagon_vertices[5]])
triangles.append([hexagon_vertices[0], hexagon_vertices[2], hexagon_vertices[5]])
triangles.append([hexagon_vertices[3], hexagon_vertices[5], hexagon_vertices[2]])

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    translated_pillar_polygon = translate_at(b, center, pillar_polygon, [], [])
    rotated_60 = rotate_at(center, 60, translated_pillar_polygon, [], [])
    rotated_120 = rotate_at(center, 120, translated_pillar_polygon, [], [])
    rotated_180 = rotate_at(center, 180, translated_pillar_polygon, [], [])
    rotated_240 = rotate_at(center, 240, translated_pillar_polygon, [], [])
    rotated_300 = rotate_at(center, 300, translated_pillar_polygon, [], [])
    pillar_polygons += [translated_pillar_polygon, rotated_60, rotated_120, rotated_180, rotated_240, rotated_300]

center = np.array([3*s + s/2, 3.0/2*s*math.sqrt(3)])
translated_vertices = translate_at(b, center, vertices, edges, edges_descriptions)
rotate_at(center, 60, translated_vertices, edges, edges_descriptions)
rotate_at(center, -60, translated_vertices, edges, edges_descriptions)

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    translated_pillar_polygon = translate_at(b, center, pillar_polygon, [], [])
    rotated_60 = rotate_at(center, 60, translated_pillar_polygon, [], [])
    rotated_300 = rotate_at(center, -60, translated_pillar_polygon, [], [])
    pillar_polygons += [translated_pillar_polygon, rotated_60, rotated_300]


center = np.array([2*s, 0.0])
translated_vertices = translate_at(b, center, vertices, edges, [])
rotate_at(center, 120, translated_vertices, edges, edges_descriptions)
rotate_at(center, 180, translated_vertices, edges, edges_descriptions)
rotate_at(center, 240, translated_vertices, edges, edges_descriptions)

# do the same with hexagon vertices
hexagon_vertices = [q1, q2]
hexagon_edges = [[0, 1]]
hexagon_edges_descriptions = []
translated_hexagon_vertices = translate_at(b, center, hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 120, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 180, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
rotate_at(center, 240, translated_hexagon_vertices, hexagon_edges, hexagon_edges_descriptions)
hexagon_vertices = hexlib.extract_vertices_from_edges_descriptions(hexagon_edges_descriptions)
triangles.append([hexagon_vertices[3], hexagon_vertices[4], hexagon_vertices[5]])
triangles.append([hexagon_vertices[3], hexagon_vertices[5], hexagon_vertices[2]])
triangles_to_be_reflected += [6, 7]

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    translated_pillar_polygon = translate_at(b, center, pillar_polygon, [], [])
    rotated_120 = rotate_at(center, 120, translated_pillar_polygon, [], [])
    rotated_180 = rotate_at(center, 180, translated_pillar_polygon, [], [])
    rotated_240 = rotate_at(center, 240, translated_pillar_polygon, [], [])
    pillar_polygons += [rotated_120, rotated_180, rotated_240]


center = np.array([3*s + s/2, s*math.sqrt(3)/2])
translated_vertices = translate_at(b, center, vertices, edges, [])
rotate_at(center, 180, translated_vertices, edges, edges_descriptions)
rotate_at(center, 240, translated_vertices, edges, edges_descriptions)
rotate_at(center, 300, translated_vertices, edges, edges_descriptions)

# same with pillar polygons
for pillar_polygon in original_pillar_polygons:
    translated_pillar_polygon = translate_at(b, center, pillar_polygon, [], [])
    rotated_180 = rotate_at(center, 180, translated_pillar_polygon, [], [])
    rotated_240 = rotate_at(center, 240, translated_pillar_polygon, [], [])
    rotated_300 = rotate_at(center, 300, translated_pillar_polygon, [], [])
    pillar_polygons += [rotated_180, rotated_240, rotated_300]

hexlib.add_new_edges(edges_descriptions, vertices, edges)

# transforming to final square shape
vertices = np.array(vertices)
transformation_matrix = 2./3 * np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
resulting_vertices = transformation_matrix * vertices.transpose()
vertices = np.asarray(resulting_vertices.transpose()) - 1.0


# deal with triangles
for index, triangle in enumerate(triangles):
    resulting_triangle = transformation_matrix * np.array(triangle).transpose()
    triangle = np.asarray(resulting_triangle.transpose()) - 1.0
    triangles[index] = triangle

# reflect triangles
reflected_triangles = []
for index in triangles_to_be_reflected:
    resulting_triangle = hexlib.reflect(135, triangles[index])
    reflected_triangles.append(resulting_triangle)

triangles += reflected_triangles

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


hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

#TODO: find, in microstructure, smaller and larger pillar tickness, so we can compute thickness void and real thickness with more precision
thickness_void = (triangle_side*pillar_area_ratio - num_pillars*thickness) / (num_pillars - 1)
min_resolution = max(2 / thickness_void, 2 / thickness)
chosen_resolution = 2 * math.pow(2, math.ceil(math.log(min_resolution) / math.log(2)))

if chosen_resolution > 1024:
        print "Resolution of " + str(chosen_resolution) + "is too big"
        print "Skipping experiment!"
        exit()

if chosen_resolution < 64:
    chosen_resolution = 64

print "Thickness void: " + str(thickness_void)
print "Minimum resolution: " + str(min_resolution)
print "Chosen resolution: " + str(chosen_resolution)

hexlib.inflate_hexagonal_box_smarter(out_wire, 0.00001, 0.0, out_mesh, incenters_thickness_pairs + pillar_nodes_custom_pairs, chosen_resolution)








