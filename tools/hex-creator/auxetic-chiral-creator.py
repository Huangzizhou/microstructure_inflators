#!/usr/bin/env python
import math
import sys
import numpy as np
import hexlib


def same_vertex(v1, v2):
    if (abs(v1[0] - v2[0]) < hexlib.Tolerance) and (v1[1] - v2[1]) < hexlib.Tolerance:
        return True

    return False


def same_edges(edge1, edge2):
    if not edge1:
        if not edge2:
            return True
        else:
            return False
    else:
        if not edge2:
            return False

    if same_vertex(edge1[0], edge2[0]) and same_vertex(edge1[1], edge2[1]):
        return True

    if same_vertex(edge1[0], edge2[1]) and same_vertex(edge1[1], edge2[0]):
        return True

    return False


def create_pillar_edge(edge):
    crossing_edge = [[-1.0, 0.0], [1.0, 0.0]]
    resulting_edge = []

    if edge[0][1] < 0.0:
        return resulting_edge

    try:
        edge_intersection_point = hexlib.edge_intersection(edge, crossing_edge)

        if edge_intersection_point[0] <= (edge[1][0] + hexlib.Tolerance):
            resulting_edge = [edge[0], edge_intersection_point]
        else:
            resulting_edge = edge

    except:
        pass

    return resulting_edge


def create_pillars_with_constant_spacing_and_thickness(line1, line2, pillar_area_ratio, num_pillars, thickness, gap):
    edges_descriptions = []
    triangles = []
    unit_vector = (line1[1] - line1[0]) / np.linalg.norm(line1[1] - line1[0])
    next_point1 = line1[0]
    next_point2 = line2[1] + (line2[0] - line2[1]) * pillar_area_ratio

    mid_line1 = (line1[0] + line1[1]) / 2

    edges_descriptions.append([line2[0], next_point2])

    for t in range(0, num_pillars):
        beginning_pillar1 = next_point1
        end_pillar1 = beginning_pillar1 + thickness * unit_vector

        beginning_pillar2 = next_point2
        end_pillar2 = beginning_pillar2 + thickness * unit_vector

        left_edge = [beginning_pillar1, end_pillar1]
        right_edge = [beginning_pillar2, end_pillar2]
        top_edge = [beginning_pillar1, beginning_pillar2]
        bottom_edge = [end_pillar1, end_pillar2]

        resulting_left_edge = create_pillar_edge(left_edge)
        resulting_right_edge = create_pillar_edge(right_edge)
        resulting_top_edge = create_pillar_edge(top_edge)
        resulting_bottom_edge = create_pillar_edge(bottom_edge)

        pillar_edges = [resulting_left_edge, resulting_right_edge, resulting_top_edge, resulting_bottom_edge]

        # create triangles
        if same_edges(top_edge, resulting_top_edge):
            if same_edges(bottom_edge, resulting_bottom_edge):
                pillar_edges.append(create_pillar_edge([beginning_pillar1, end_pillar2]))

                triangles.append([beginning_pillar1, end_pillar2, end_pillar1])
                triangles.append([beginning_pillar1, end_pillar2, beginning_pillar2])

            else:
                pillar_edges.append([beginning_pillar1, resulting_bottom_edge[1]])
                pillar_edges.append([resulting_bottom_edge[1], beginning_pillar2])
                pillar_edges.append([resulting_bottom_edge[1], [beginning_pillar2[0], 0.0]])

                triangles.append([beginning_pillar1, end_pillar1, resulting_bottom_edge[1]])
                triangles.append([beginning_pillar1, resulting_bottom_edge[1], beginning_pillar2])
                triangles.append([beginning_pillar2, resulting_bottom_edge[1], [beginning_pillar2[0], 0.0]])
        else:
            if not resulting_top_edge:
                pass
            else:
                if not resulting_bottom_edge:
                    pillar_edges.append([resulting_top_edge[1], [beginning_pillar1[0], 0.0]])

                    triangles.append([resulting_top_edge[1], [beginning_pillar1[0], 0.0], beginning_pillar1])
                else:
                    pillar_edges.append([beginning_pillar1, resulting_bottom_edge[1]])
                    pillar_edges.append([resulting_top_edge[1], resulting_bottom_edge[1]])

                    triangles.append([beginning_pillar1, resulting_bottom_edge[1], resulting_top_edge[1]])
                    triangles.append([beginning_pillar1, end_pillar1, resulting_bottom_edge[1]])

        # edges_descriptions.append(create_pillar_edge([beginning_pillar1, end_pillar1]))
        # edges_descriptions.append(create_pillar_edge([beginning_pillar2, end_pillar2]))
        # edges_descriptions.append(create_pillar_edge([beginning_pillar1, beginning_pillar2]))
        # edges_descriptions.append(create_pillar_edge([end_pillar1, end_pillar2]))
        # edges_descriptions.append(create_pillar_edge([beginning_pillar1, end_pillar2])) # create triangles

        next_point1 = end_pillar1 + gap * unit_vector
        next_point2 = end_pillar2 + gap * unit_vector

        pillar_edges.append(create_pillar_edge([end_pillar1, next_point1]))
        pillar_edges.append(create_pillar_edge([end_pillar2, next_point2]))

        for pillar_edge in pillar_edges:
            if not pillar_edge:
                pass
            else:
                edges_descriptions.append(pillar_edge)

    # TODO: only add if different point. Don't want edge with zero length

    final_edge = create_pillar_edge([end_pillar1, mid_line1])
    if not final_edge:
        pass
    else:
        edges_descriptions.append(final_edge)

    return edges_descriptions, triangles


def subdivide_segment_with_constant_pillar_spacing_and_thickness(a, b, num_pillars, thickness, gap):
    edges_descriptions = []
    unit_vector = (b - a) / np.linalg.norm(b - a)
    next_point = a

    for t in range(0, num_pillars):
        beginning_pillar = next_point
        middle_pillar = beginning_pillar + thickness / 2 * unit_vector
        end_pillar = beginning_pillar + thickness * unit_vector

        edges_descriptions.append([beginning_pillar, middle_pillar])
        edges_descriptions.append([middle_pillar, end_pillar])
        edges_descriptions.append([beginning_pillar, end_pillar])

        if t < (num_pillars - 1):
            next_point = end_pillar + gap * unit_vector
            edges_descriptions.append([end_pillar, next_point])

    return edges_descriptions


def simplex_to_whole_parallelogram(simplex_vertices, simplex_edges, parallelogram_side):
    l = parallelogram_side / 2
    vertices = list(simplex_vertices)
    edges = list(simplex_edges)

    initial_vertices = list(vertices)
    initial_edges = list(edges)

    # Operation 1a: create new vertices by rotating about b 120 degrees
    b = np.array([l, l * math.sqrt(3) / 3.0])
    translated_vertices = vertices - b
    rotated_vertices = hexlib.rotate(120, translated_vertices)

    rotated_edges = list(edges)
    rotated_vertices += b
    hexlib.add_new_vertices_and_edges(rotated_vertices, rotated_edges, vertices, edges)

    # Operation 1b: create new vertices by rotating initial structure about b -120 degrees
    translated_vertices = initial_vertices - b
    rotated_vertices = hexlib.rotate(240, translated_vertices)

    rotated_edges = list(initial_edges)
    rotated_vertices += b
    hexlib.add_new_vertices_and_edges(rotated_vertices, rotated_edges, vertices, edges)

    # Operation 2: create new vertices by rotating current structure about the top vertex 60 degrees
    top_vertex = np.array([l, l * math.sqrt(3)])
    translated_vertices = vertices - top_vertex
    rotated_vertices = hexlib.rotate(60, translated_vertices)

    rotated_edges = list(edges)
    rotated_vertices += top_vertex
    hexlib.add_new_vertices_and_edges(rotated_vertices, rotated_edges, vertices, edges)

    # Operation 5: transform parallelogram to square
    vertices = np.array(vertices)
    transformation_matrix = np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
    resulting_vertices = transformation_matrix * vertices.transpose()
    vertices = np.asarray(resulting_vertices.transpose())

    # Finally, transpose to origin
    vertices -= l

    return vertices, edges


def simplex_polygon_to_whole_parallelogram(simplex_polygon, parallelogram_side):
    l = parallelogram_side / 2
    vertices = list(simplex_polygon)
    polygons = [simplex_polygon]

    initial_vertices = list(vertices)

    # Operation 1a: create new vertices by rotating about b 120 degrees
    b = np.array([l, l * math.sqrt(3) / 3.0])
    translated_vertices = vertices - b
    rotated_vertices = hexlib.rotate(120, translated_vertices)
    rotated_vertices += b

    polygons.append(rotated_vertices)  # polygons now has 2 pols

    # Operation 1b: create new vertices by rotating initial structure about b -120 degrees
    translated_vertices = initial_vertices - b
    rotated_vertices = hexlib.rotate(240, translated_vertices)
    rotated_vertices += b

    polygons.append(rotated_vertices)  # polygons now has 3 pols

    # Operation 2: create new vertices by rotating current structure about the top vertex 60 degrees
    polygons_copy = list(polygons)
    for index, pol in enumerate(polygons_copy):
        vertices = list(pol)
        top_vertex = np.array([l, l * math.sqrt(3)])
        translated_vertices = vertices - top_vertex
        rotated_vertices = hexlib.rotate(60, translated_vertices)
        rotated_vertices += top_vertex

        polygons.append(rotated_vertices)

    # Operation 5: transform parallelogram to square
    transformation_matrix = np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
    for index, pol in enumerate(polygons):
        vertices = np.array(pol)
        resulting_vertices = transformation_matrix * vertices.transpose()
        polygons[index] = np.asarray(resulting_vertices.transpose() - l)

    return polygons


if len(sys.argv) != 7:
    print "usage: ./auxetic-chiral-creator.py <size of triangles> <num of pillars> <relative size of pillar area (chirality)> <relative thickness (of pillars)> " \
          "<output wire> <output mesh>"
    print "example: ./auxetic-chiral-creator.py 0.8 3 0.6 0.9 output.wire output.msh"
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
parallelogram_side = 2.0
l = parallelogram_side / 2.0

# graph to be constructed
vertices = []
edges = []

# define important vertices of simplex used to build entire parallelogram structure
origin = np.array([0, 0])
a = origin.copy()
b = np.array([l, l * math.sqrt(3) / 3.0])
b_yneg = np.array([l, -l * math.sqrt(3) / 3.0])
c = np.array([2 * l, 0])
d = np.array([l, 0])

# define vertices on the left side of the simplex based on parameters
height = l * math.sqrt(3) / 3.0
s = 2 * height
triangle_side = triangle_side_ratio * s

q1 = triangle_side_ratio * (b - a)
q2 = triangle_side_ratio * (d - a)
q3 = triangle_side_ratio * (b_yneg - a)

# pillar area: extremes are q1 and q4 on left side
q4 = q1 + (q3 - q1) * pillar_area_ratio

# individual pillar area:
thickness = thickness_ratio * (pillar_area_ratio * s / num_pillars)
spacing = (pillar_area_ratio * triangle_side_ratio * s - thickness) / (num_pillars - 1) - thickness

# compute points on the left side
# edges_descriptions_left = subdivide_segment_with_constant_pillar_spacing_and_thickness(q1, q4, num_pillars, thickness, spacing)
# connecting_vertices_left = hexlib.extract_vertices_from_edges_descriptions(edges_descriptions_left)

# define vertices on the right side of the simplex based on parameters
w1 = triangle_side_ratio * (b - c) + c
w2 = triangle_side_ratio * (d - c) + c
w3 = triangle_side_ratio * (b_yneg - c) + c

# pillar area: extremes are q1 and q4 on right side
# w4 = w3 + (w1 - w3) * pillar_area_ratio

print "Constructing " + out_wire + " ..."

edge_descriptions, triangles = create_pillars_with_constant_spacing_and_thickness([q1, q3], [w1, w3], pillar_area_ratio,
                                                                                  num_pillars, thickness, spacing)

# add edges to fill big inner hexagons
edge_descriptions.append([a, q1])
edge_descriptions.append([a, q2])
edge_descriptions.append([c, w1])
edge_descriptions.append([c, w2])

triangles.append([a, q1, q2])
triangles.append([c, w1, w2])

hexlib.add_new_edges(edge_descriptions, vertices, edges)

vertices, edges = simplex_to_whole_parallelogram(vertices, edges, parallelogram_side)

vertices = vertices.tolist()

parallelogram_triangles = []
for triangle in triangles:
    parallelogram_triangles += simplex_polygon_to_whole_parallelogram(triangle, parallelogram_side)

triangles_incenters_thickness_pairs = hexlib.add_polygons_incenters(parallelogram_triangles, vertices, edges)

hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

hexlib.inflate_hexagonal_box_smarter(out_wire, 0.00001, 0.00001, out_mesh, triangles_incenters_thickness_pairs)