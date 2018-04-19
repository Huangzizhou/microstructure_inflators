#!/usr/bin/env python
import math
import sys
import numpy as np
import hexlib

if len(sys.argv) != 6:
    print "usage: ./auxetic-star-creator.py <ab positioning> <ac positioning> <thickeness> " \
          "<output wire> <output mesh>"
    print "example: ./auxetic-star-creator.py 0.4 0.7 0.5 output.wire output.msh"
    sys.exit(-1)

#          . b
#     . '  |
# a._______| c
#
ab_param = float(sys.argv[1])
ac_param = float(sys.argv[2])
thickness = float(sys.argv[3])
out_wire = sys.argv[4]
out_mesh = sys.argv[5]
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
vertices.append(c)

# fundamental edges
edges.append([0, 1])
edges.append([1, 2])

# create topology in parallelogram
vertices, edges = hexlib.simplex_to_whole_square(vertices, edges)

# print wire output
hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

hexlib.inflate_hexagonal_box_smarter(out_wire, thickness, 0.001, out_mesh)
