#!/usr/bin/env python
import math
import sys
import numpy as np
import hexlib


if len(sys.argv) != 7:
    print "usage: ./hexa-many-pillars-russian-parameters.py <size of triangles> <num of pillars> <relative size of pillar area (chirality)> <relative thickness (of pillars)> " \
          "<output wire> <output mesh>"
    print "example: ./hexa-many-pillars.py 0.8 3 0.6 0.9 output.wire output.msh"
    sys.exit(-1)

# . b
#     . '  |
# a._______| c
#       d
triangle_side_ratio = float(sys.argv[1])
num_pillars = int(sys.argv[2])
pillar_area_ratio = float(sys.argv[3])
thickness_ratio = float(sys.argv[4])
out_wire = sys.argv[5]
out_mesh = sys.argv[6]
parallelogram_side = 2.0
l = parallelogram_side / 2.0

s = 2 / math.sqrt(3)

triangle_side = triangle_side_ratio * s * math.sqrt(3)
thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars)

num_upper_intervals = num_pillars - 1
if num_upper_intervals > 0:
    upper_spacing = (triangle_side - thickness) / (num_pillars-1)
    num_upper_points = int(math.floor(num_pillars / 2.0) + 1)

# graph to be constructed
vertices = []
edges = []

# define important vertices of simplex used to build entire parallelogram structure
origin = np.array([0, 0])
a = origin.copy()
b = np.array([l, l * math.sqrt(3) / 3.0])
c = np.array([l, 0])

# positions can be inferred from the triangle side chosen
radius = triangle_side * math.sqrt(3) / 6
f = b - [0, radius]

# define vertices on simplex based on parameters
p1 = f - [triangle_side / 2.0, 0]
p2 = f

upper_triangles = hexlib.simplex_polygon_to_whole_parallelogram([b, p1, p2], parallelogram_side)

print "Constructing " + out_wire + " ..."

# Add vertices and edges of the first simplex
vertices.append(b)
vertices.append(p1)
p1_offset = [p1[0] + thickness/2, p1[1]]
vertices.append(p1_offset)

if (num_pillars % 2) == 1:
    vertices.append(c)

# fundamental edges
edges.append([0, 1])
edges.append([1, 2])

# creating vertices and edges between p1 and p2
edges_descriptions = []

if num_upper_intervals > 0:
    hexlib.create_subdivided_segment_with_constant_spacing(p1_offset, p2, num_upper_points, upper_spacing, edges_descriptions)
    vertices_from_p1_to_p2 = hexlib.extract_vertices_from_edges_descriptions(edges_descriptions)
    hexlib.add_new_edges(edges_descriptions, vertices, edges)
else:
    vertices_from_p1_to_p2 = [p2]
    vertices.append(p2)

# adding edge between p2, b and p2,c
edges.append([len(vertices)-1, 0])
if (num_pillars % 2) == 1:
    edges.append([len(vertices)-1, 2])


if num_upper_intervals > -1:
    pillar_triangles = []
    for index, top_point in enumerate(vertices_from_p1_to_p2):
        if index == (len(vertices_from_p1_to_p2)-1):
            if (num_pillars % 2) == 0: # if even number of nodes, last one does not have pillars
                continue
            else:
                a1 = np.array([top_point[0] - thickness / 2, top_point[1]])
                a2 = np.array([top_point[0], top_point[1]])

                b1 = np.array([top_point[0] - thickness / 2, 0])
                b2 = np.array([top_point[0], 0])
        else:
            a1 = np.array([top_point[0] - thickness/2, top_point[1]])
            a2 = np.array([top_point[0] + thickness / 2, top_point[1]])

            b1 = np.array([top_point[0] - thickness/2, 0])
            b2 = np.array([top_point[0] + thickness/2, 0])

        triangle_left = [a1, b1, b2]
        triangle_right = [a1, b2, a2]

        edges_descriptions += hexlib.polygon_to_edges_descriptions(triangle_left)
        edges_descriptions += hexlib.polygon_to_edges_descriptions(triangle_right)
        hexlib.add_new_edges(edges_descriptions, vertices, edges)

        pillar_triangles += hexlib.simplex_polygon_to_whole_parallelogram(triangle_left, parallelogram_side)
        pillar_triangles += hexlib.simplex_polygon_to_whole_parallelogram(triangle_right, parallelogram_side)
else:
    pillar_triangles = []

# create topology in parallelogram
vertices, edges = hexlib.simplex_to_whole_parallelogram(vertices, edges, parallelogram_side)
vertices = vertices.tolist()

upper_incenters_thickness_pairs = hexlib.add_polygons_incenters(upper_triangles, vertices, edges)
pillar_triangles_incenters_thickness_pairs = hexlib.add_polygons_incenters(pillar_triangles, vertices, edges)

# print wire output
hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

#hexlib.inflate_hexagonal_box_smarter(out_wire, 0.005, 0.001, out_mesh, upper_incenters_thickness_pairs + pillar_triangles_incenters_thickness_pairs)
hexlib.inflate_hexagonal_box_smarter(out_wire, 0.00001, 0.00001, out_mesh, upper_incenters_thickness_pairs + pillar_triangles_incenters_thickness_pairs)
