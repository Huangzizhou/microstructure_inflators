#!/usr/bin/env python
import math
import sys

import numpy as np

import hexlib
from hexlib import polygon_centroid, min_distance_point_line, simplex_vertices_to_whole_parallelogram, \
    simplex_to_whole_parallelogram

if len(sys.argv) != 6:
    print "usage: ./six-leaf-clover-creator.py <ab positioning> <bd positioning> <cb positioning> " \
          "<output wire> <output mesh>"
    print "example: ./six-leaf-clover-creator.py 0.8 0.6 0.7 output.wire output.msh"
    sys.exit(-1)

# . b
#     . '  |
# a.____.__| c
#       d
ab_param = float(sys.argv[1])
bd_param = float(sys.argv[2])
cb_param = float(sys.argv[3])
out_wire = sys.argv[4]
out_mesh = sys.argv[5]
parallelogram_side = 2.0
l = parallelogram_side / 2.0

# graph to be constructed
vertices = []
edges = []

# define important vertices of simplex used to build entire parallelogram structure
origin = np.array([0, 0])
a = origin.copy()
b = np.array([l, l * math.sqrt(3) / 3.0])
c = np.array([l, 0])
d = np.array([2.0 * l / 3, 0])

# define vertices on simplex based on parameters
p1 = (1.0 - ab_param) * b + ab_param * a
p2 = (1.0 - bd_param) * d + bd_param * b
p3 = (1.0 - cb_param) * b + cb_param * c

# auxiliary vertices to fill space
midpoint_left = polygon_centroid([p1, d, p2])
midpoint_right = polygon_centroid([p2, d, p3])
p1_middle = [p1[0], p1[1]/2]
p3_middle = [p3[0], p3[1]/2]

print "Constructing " + out_wire + " ..."

# Add vertices and edges of the first simplex
vertices.append(p1)
vertices.append(p2)
vertices.append(p3)
vertices.append(a)
vertices.append(c)
vertices.append(d)

vertices.append(midpoint_left)
vertices.append(p1_middle)
vertices.append(p3_middle)
vertices.append(midpoint_right)

# fundamental edges
edges.append([0, 1])
edges.append([1, 2])
edges.append([3, 5])
edges.append([4, 5])
edges.append([1, 5])

# edges for filling space
edges.append([0, 3])
edges.append([3, 7])
edges.append([7, 6])
edges.append([6, 8])
edges.append([6, 1])
edges.append([7, 8])
edges.append([9, 2])
edges.append([9, 1])

# create edges and vertices to fill in the polygons
midpoint_left_vertices = simplex_vertices_to_whole_parallelogram([midpoint_left], parallelogram_side)
midpoint_right_vertices = simplex_vertices_to_whole_parallelogram([midpoint_right], parallelogram_side)
p1_middle_vertices = simplex_vertices_to_whole_parallelogram([p1_middle], parallelogram_side)
p3_middle_vertices = simplex_vertices_to_whole_parallelogram([p3_middle], parallelogram_side)
d_vertices = simplex_vertices_to_whole_parallelogram([d], parallelogram_side)

# create topology in parallelogram
vertices, edges = simplex_to_whole_parallelogram(vertices, edges, parallelogram_side)

# print wire output
hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

thickness_midpoint_left = 1.2 * min_distance_point_line(midpoint_left, [p1, p2])
thickness_midpoint_right = 1.15 * min_distance_point_line(midpoint_right, [p2, p3])
thickness_p1_middle = p1[1] * 0.45
thickness_p3_middle = p3[1] * 0.45

hexlib.inflate_hexagonal_box_smarter(out_wire, 0.02, 0.001,  out_mesh, [midpoint_left_vertices, thickness_midpoint_left],
                                     [midpoint_right_vertices, thickness_midpoint_right], [p1_middle_vertices, thickness_p1_middle],
                                     [p3_middle_vertices, thickness_p3_middle], [d_vertices, 0.2])