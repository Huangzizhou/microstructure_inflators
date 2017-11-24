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


if len(sys.argv) != 5:
    print "usage: ./hexa-diamond-creator.py <size of triangles> <thickeness of diamond pillar> " \
          "<output wire> <output mesh>"
    print "example: ./hexa-diamond-creator.py 1.0 0.3 output.wire output.msh"
    sys.exit(-1)

# . b
#     . '  |
# a._______| c
#       d
triangle_side = float(sys.argv[1])
pillars_thickness = float(sys.argv[2])
out_wire = sys.argv[3]
out_mesh = sys.argv[4]
parallelogram_side = 2.0
l = parallelogram_side / 2.0
s = triangle_side

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
p3 = c - [pillars_thickness, 0]

upper_triangles = hexlib.simplex_polygon_to_whole_parallelogram([b, p1, p2], parallelogram_side)
lower_triangles = hexlib.simplex_polygon_to_whole_parallelogram([p2, p3, c], parallelogram_side)

print "Constructing " + out_wire + " ..."

# Add vertices and edges of the first simplex
vertices.append(b)
vertices.append(p1)
vertices.append(p2)
vertices.append(p3)
vertices.append(c)

# fundamental edges
edges.append([0, 1])
edges.append([1, 2])
edges.append([2, 3])
edges.append([3, 4])
edges.append([4, 2])
edges.append([2, 0])

# create edges and vertices to fill in the polygons
b_vertices = hexlib.simplex_vertices_to_whole_parallelogram([b], parallelogram_side)

# create topology in parallelogram
vertices, edges = hexlib.simplex_to_whole_parallelogram(vertices, edges, parallelogram_side)
vertices = vertices.tolist()

upper_incenters_thickness_pairs = add_polygons_incenters(upper_triangles, vertices, edges)
lower_incenters_thickness_pairs = add_polygons_incenters(lower_triangles, vertices, edges)

# print wire output
hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

centroid_thickness = hexlib.min_distance_point_line()

hexlib.inflate_hexagonal_box_smarter(out_wire, 0.005, 0.001, out_mesh, upper_incenters_thickness_pairs + lower_incenters_thickness_pairs)
