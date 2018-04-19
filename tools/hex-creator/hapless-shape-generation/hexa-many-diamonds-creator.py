#!/usr/bin/env python
import math
import sys
import numpy as np
import hexlib


def add_polygons_incenters(polygons, vertices, edges):
    incenters_thickness_pairs = []

    for index, pol in enumerate(polygons):
        incenter = hexlib.triangle_incenter(pol)
        thickness = min(hexlib.min_distance_point_line(incenter, [pol[0], pol[1]]),
                        hexlib.min_distance_point_line(incenter, [pol[1], pol[2]]),
                        hexlib.min_distance_point_line(incenter, [pol[2], pol[0]]))

        incenters_thickness_pairs.append([[incenter], thickness])

        edges_descriptions = []
        for i in range(0, len(pol)):
            a = pol[i].tolist()
            b = incenter
            edges_descriptions.append([a, b])

        hexlib.add_new_edges(edges_descriptions, vertices, edges)

    return incenters_thickness_pairs


if len(sys.argv) != 6:
    print "usage: ./hexa-diamond-creator.py <size of triangles> <num of pillars> <thickeness of diamond pillar> " \
          "<output wire> <output mesh>"
    print "example: ./hexa-diamond-creator.py 1.0 3 0.3 output.wire output.msh"
    sys.exit(-1)

# . b
#     . '  |
# a._______| c
#       d
triangle_side = float(sys.argv[1])
num_pillars = int(sys.argv[2])
pillars_thickness = float(sys.argv[3])
out_wire = sys.argv[4]
out_mesh = sys.argv[5]
parallelogram_side = 2.0
l = parallelogram_side / 2.0
s = triangle_side

num_upper_intervals = num_pillars - 1
upper_spacing = (s - pillars_thickness) / (num_pillars-1)
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
radius = s * math.sqrt(3) / 6
f = b - [0, radius]

# define vertices on simplex based on parameters
p1 = f - [s / 2.0, 0]
p2 = f


upper_triangles = hexlib.simplex_polygon_to_whole_parallelogram([b, p1, p2], parallelogram_side)

print "Constructing " + out_wire + " ..."

# Add vertices and edges of the first simplex
vertices.append(b)
vertices.append(p1)
p1_offset = [p1[0] + pillars_thickness/2, p1[1]]
vertices.append(p1_offset)

if (num_pillars % 2) == 1:
    vertices.append(c)

# fundamental edges
edges.append([0, 1])
edges.append([1, 2])

# creating vertices and edges between p1 and p2
edges_descriptions = []
hexlib.create_subdivided_segment_with_constant_spacing(p1_offset, p2, num_upper_points, upper_spacing, edges_descriptions)
vertices_from_p1_to_p2 = hexlib.extract_vertices_from_edges_descriptions(edges_descriptions)
hexlib.add_new_edges(edges_descriptions, vertices, edges)

# adding edge between p2, b and p2,c
edges.append([len(vertices)-1, 0])
if (num_pillars % 2) == 1:
    edges.append([len(vertices)-1, 2])

pillar_triangles = []
for index, top_point in enumerate(vertices_from_p1_to_p2):
    if index == (len(vertices_from_p1_to_p2)-1) and (num_pillars % 2) == 0: # if even number of nodes, last one does not have pillars
        continue

    b1 = np.array([top_point[0] - pillars_thickness/2, 0])
    b2 = np.array([top_point[0] + pillars_thickness/2, 0])

    triangle = [top_point, b1, b2]
    edges_descriptions = [[top_point, b1], [top_point, b2], [b1, b2]]
    hexlib.add_new_edges(edges_descriptions, vertices, edges)

    pillar_triangles += hexlib.simplex_polygon_to_whole_parallelogram(triangle, parallelogram_side)


# create topology in parallelogram
vertices, edges = hexlib.simplex_to_whole_parallelogram(vertices, edges, parallelogram_side)
vertices = vertices.tolist()

upper_incenters_thickness_pairs = add_polygons_incenters(upper_triangles, vertices, edges)
pillar_triangles_incenters_thickness_pairs = add_polygons_incenters(pillar_triangles, vertices, edges)

# print wire output
hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

hexlib.inflate_hexagonal_box_smarter(out_wire, 0.005, 0.001, out_mesh, upper_incenters_thickness_pairs + pillar_triangles_incenters_thickness_pairs)
