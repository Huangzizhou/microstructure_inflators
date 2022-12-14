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
    print "usage: ./hexa-many-pillars.py <size of triangles> <num of pillars> <thickeness of pillars> " \
          "<output wire> <output mesh>"
    print "example: ./hexa-many-pillars.py 1.0 3 0.3 output.wire output.msh"
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
if num_upper_intervals > 0:
    upper_spacing = (s - pillars_thickness) / (num_pillars - 1)
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
p1_offset = [p1[0] + pillars_thickness / 2, p1[1]]
vertices.append(p1_offset)

if (num_pillars % 2) == 1:
    vertices.append(c)

# fundamental edges
edges.append([0, 1])
edges.append([1, 2])

# creating vertices and edges between p1 and p2
edges_descriptions = []

if num_upper_intervals > 0:
    hexlib.create_subdivided_segment_with_constant_spacing(p1_offset, p2, num_upper_points, upper_spacing,
                                                           edges_descriptions)
    vertices_from_p1_to_p2 = hexlib.extract_vertices_from_edges_descriptions(edges_descriptions)
    hexlib.add_new_edges(edges_descriptions, vertices, edges)
else:
    vertices_from_p1_to_p2 = [p2]
    vertices.append(p2)

# adding edge between p2, b and p2,c
edges.append([len(vertices) - 1, 0])
if (num_pillars % 2) == 1:
    edges.append([len(vertices) - 1, 2])

if num_upper_intervals > -1:
    pillar_triangles = []
    for index, top_point in enumerate(vertices_from_p1_to_p2):
        if index == (len(vertices_from_p1_to_p2) - 1):
            if (num_pillars % 2) == 0:  # if even number of nodes, last one does not have pillars
                continue
            else:
                a1 = np.array([top_point[0] - pillars_thickness / 2, top_point[1]])
                a2 = np.array([top_point[0], top_point[1]])

                b1 = np.array([top_point[0] - pillars_thickness / 2, 0])
                b2 = np.array([top_point[0], 0])
        else:
            a1 = np.array([top_point[0] - pillars_thickness / 2, top_point[1]])
            a2 = np.array([top_point[0] + pillars_thickness / 2, top_point[1]])

            b1 = np.array([top_point[0] - pillars_thickness / 2, 0])
            b2 = np.array([top_point[0] + pillars_thickness / 2, 0])

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

upper_incenters_thickness_pairs = add_polygons_incenters(upper_triangles, vertices, edges)
pillar_triangles_incenters_thickness_pairs = add_polygons_incenters(pillar_triangles, vertices, edges)

# print wire output
hexlib.create_wire(vertices, edges, out_wire)

print "Inflating ..."

# Computing void thickness and necessary resolution
thickness_void = (triangle_side - num_pillars * pillars_thickness) / (num_pillars - 1)
min_resolution = max(2 / thickness_void, 2 / pillars_thickness)
chosen_resolution = math.pow(2, math.ceil(math.log(min_resolution) / math.log(2)))

if chosen_resolution > 1024:
        print "Resolution of " + str(chosen_resolution) + "is too big"
        print "Skipping experiment!"
        exit()

if chosen_resolution < 64:
    chosen_resolution = 64

print "Thickness void: " + str(thickness_void)
print "Minimum resolution: " + str(min_resolution)
print "Chosen resolution: " + str(chosen_resolution)

hexlib.inflate_hexagonal_box_smarter(out_wire, 0.00001, 0.00001, out_mesh,
                                     upper_incenters_thickness_pairs + pillar_triangles_incenters_thickness_pairs,
                                     chosen_resolution)
