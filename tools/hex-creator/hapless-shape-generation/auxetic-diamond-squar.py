#!/usr/bin/env python
import math
import sys
import numpy as np
import hexlib

DIAMOND_TILING = False

if len(sys.argv) != 7:
    print "usage: ./auxetic-diamond-squar-creator.py <ab positioning> <ac positioning> <star thickeness> <tiling thickness>" \
          "<output wire> <output mesh>"
    print "example: ./auxetic-diamond-squar-creator.py 0.4 0.7 0.5 output.wire output.msh"
    sys.exit(-1)

#          . b
#     . '  |
# a._______| c
#
ab_param = float(sys.argv[1])
ac_param = float(sys.argv[2])
thickness = float(sys.argv[3])
tiling_thickness = float(sys.argv[4])
out_wire = sys.argv[5]
out_mesh = sys.argv[6]
square_side = 2.0
l = square_side / 2.0

# graph to be constructed
vertices = []
edges = []

# define important vertices of simplex used to build entire parallelogram structure
origin = np.array([0, 0])
a = origin.copy()
b = np.array([l, l])
c = np.array([l, 0])
d = np.array([2.0 * l / 3, 0])

# define vertices on simplex based on parameters
p1 = (1.0 - ab_param) * b + ab_param * a
p2 = (1.0 - ac_param) * c + ac_param * a

print "Constructing " + out_wire + " ..."

# Add vertices and edges of the first simplex
vertices.append(p1)
vertices.append(p2)
#vertices.append(c)

# fundamental edges
#edges.append([1, 2])


# top diamond
diamond = hexlib.create_diamond(p1, p2, thickness)
edges_descriptions = hexlib.polygon_to_edges_descriptions(diamond)

upper_diamond_centroid = hexlib.polygon_centroid(diamond)

# add centroid to edges
edges_descriptions += [[upper_diamond_centroid, diamond[0]], [upper_diamond_centroid, diamond[1]], [upper_diamond_centroid, diamond[2]], [upper_diamond_centroid, diamond[3]]]

# add new edges
hexlib.add_new_edges(edges_descriptions, vertices, edges)


if DIAMOND_TILING:
    # lover diamond (only part of it)
    lower_diamond_top = np.array([(p2[0] + c[0]) / 2, tiling_thickness / 2.0])
    lower_triangle = [lower_diamond_top, p2, c]
    edges_descriptions = hexlib.polygon_to_edges_descriptions(lower_triangle)

    lower_diamond_centroid = [(p2[0] + c[0])/2.0, 0.0]

    # add centroid to edges
    edges_descriptions += [[lower_diamond_centroid, lower_triangle[0]], [lower_diamond_centroid, lower_triangle[1]], [lower_diamond_centroid, lower_triangle[2]]]

    # add new edges
    hexlib.add_new_edges(edges_descriptions, vertices, edges)

    upper_diamond_centroids = hexlib.simplex_vertices_to_whole_square([upper_diamond_centroid])
    lower_diamond_centroids = hexlib.simplex_vertices_to_whole_square([lower_diamond_centroid])
    upper_centroids_thickness = hexlib.min_distance_point_line(upper_diamond_centroid, [diamond[0], diamond[1]])
    lower_centroids_thickness = hexlib.min_distance_point_line(lower_diamond_centroid, [p2, lower_diamond_top])

    custom_vertices = [[upper_diamond_centroids, upper_centroids_thickness], [lower_diamond_centroids, lower_centroids_thickness]]
else:
    vertices.append(c)
    edges.append([1, len(vertices) - 1])

    lower_bar_vertices = hexlib.simplex_vertices_to_whole_square([p2, c])
    lower_bar_thickness = tiling_thickness / 2.0

    diamond_centroids = hexlib.simplex_vertices_to_whole_square([upper_diamond_centroid])
    centroids_thickness = hexlib.min_distance_point_line(upper_diamond_centroid, [diamond[0], diamond[1]])

    custom_vertices = [[diamond_centroids, centroids_thickness], [lower_bar_vertices, lower_bar_thickness]]



# create topology in square
vertices, edges = hexlib.simplex_to_whole_square(vertices, edges)

# print wire output
hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

hexlib.inflate_hexagonal_box_smarter(out_wire, 0.005, 0.001, out_mesh, custom_vertices)
