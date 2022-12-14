#!/usr/bin/env python
import math
import os
import sys

import numpy as np

import hexlib
from hexlib import inflate_hexagonal_box_smarter

if len(sys.argv) != 7:
    print("usage: ./hex-creator.py <size of triangles> <num of pillars> <relative size of pillar area (chirality)> <relative thickness (of pillars)> " \
          "<output wire> <output mesh>")
    print("example: ./hex-creator.py 0.8 3 0.6 0.9 output.wire output.msh")
    sys.exit(-1)

# set scripts directory, so it can find all necessary files:
pathname = os.path.dirname(sys.argv[0])
hexlib.script_directory = os.path.abspath(pathname)

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

#s = 2 / math.sqrt(3)

triangle_side = triangle_side_ratio * 2
thickness = thickness_ratio * (pillar_area_ratio * triangle_side / num_pillars) / 2 # use the radius
# graph constructed on left triangle
vertices = []
edges = []

# define important vertices of parallelogram
origin = np.array([0, 0])
triangle_centroid = np.array([l, l * math.sqrt(3) / 3.0])

print("Constructing " + out_wire + " ...")

# First side: create pillars of the left side (counter-clock orientation)
#     /.
#    /  .
#   /.....
parallelogram_start = np.array([(l + triangle_side / 2.0) / 2.0, (l + triangle_side / 2.0) * math.sqrt(3) / 2.0])
parallelogram_end = np.array([(l - triangle_side / 2.0) / 2.0, (l - triangle_side / 2.0) * math.sqrt(3) / 2.0])

radius = triangle_side * math.sqrt(3) / 6
f = triangle_centroid - [0, radius]
triangle_start = f + [0, math.sqrt(3) * triangle_side / 2.0]
triangle_end = f - [triangle_side / 2.0, 0]

unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

parallelogram_start_offset = parallelogram_start + thickness * unit_vector
parallelogram_end_offset = parallelogram_end - thickness * unit_vector
triangle_start_offset = triangle_start + thickness * unit_vector
triangle_end_offset = triangle_end - thickness * unit_vector

# add new vertices and edges to current sets
new_edges_description = hexlib.create_pillars([parallelogram_start_offset, parallelogram_end_offset],
                                              [triangle_start_offset, triangle_end_offset], num_pillars)

# save pillar nodes
pillar_nodes = []
for edge_description in new_edges_description:
    pillar_nodes.append(edge_description[0])
    pillar_nodes.append(edge_description[1])

#new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
hexlib.add_new_edges(new_edges_description, vertices, edges)

# Second side: create pillars of the bottom side (counter-clock orientation)
#     ..
#    .  .
#   ______
parallelogram_start = np.array([l - triangle_side / 2.0, 0])
parallelogram_end = np.array([l + triangle_side / 2.0, 0])

triangle_start = f - [triangle_side / 2.0, 0]
triangle_end = f + [triangle_side / 2.0, 0]

pillar_example = [[parallelogram_start, triangle_start]]

unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

parallelogram_start_offset = parallelogram_start + thickness * unit_vector
parallelogram_end_offset = parallelogram_end - thickness * unit_vector
triangle_start_offset = triangle_start + thickness * unit_vector
triangle_end_offset = triangle_end - thickness * unit_vector

# add new vertices and edges to current sets
new_edges_description = hexlib.create_pillars([parallelogram_start_offset, parallelogram_end_offset],
                                              [triangle_start_offset, triangle_end_offset],
                                              num_pillars)

# save pillar nodes
for edge_description in new_edges_description:
    pillar_nodes.append(edge_description[0])
    pillar_nodes.append(edge_description[1])

#new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
hexlib.add_new_edges(new_edges_description, vertices, edges)

# Third side: create pillars of the right side (counter-clock orientation)
#     .\
#    .  \
#   .....\
parallelogram_start = np.array([2 * l - (l - triangle_side / 2.0) / 2, (l - triangle_side / 2.0) * math.sqrt(3) / 2])
parallelogram_end = np.array([2 * l - (l + triangle_side / 2.0) / 2, (l + triangle_side / 2.0) * math.sqrt(3) / 2])

triangle_start = f + [triangle_side / 2.0, 0]
triangle_end = f + [0, math.sqrt(3) * triangle_side / 2.0]

unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

parallelogram_start_offset = parallelogram_start + thickness * unit_vector
parallelogram_end_offset = parallelogram_end - thickness * unit_vector
triangle_start_offset = triangle_start + thickness * unit_vector
triangle_end_offset = triangle_end - thickness * unit_vector

# add new vertices and edges to current sets
new_edges_description = hexlib.create_pillars([parallelogram_start_offset, parallelogram_end_offset],
                                              [triangle_start_offset, triangle_end_offset],
                                              num_pillars)

# save nodes on hypotenuse
hypotenuse_nodes = []
for edge_description in new_edges_description:
    hypotenuse_nodes.append(edge_description[0])
    hypotenuse_nodes.append(edge_description[1])

#new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
hexlib.add_new_edges(new_edges_description, vertices, edges)


# Now, adding edges on the triangle
triangle = [f + [0, math.sqrt(3) * triangle_side / 2.0], f - [triangle_side / 2.0, 0], f + [triangle_side / 2.0, 0]]
new_edges_description = hexlib.create_triangle_edges(triangle, num_pillars, thickness, thickness)
#new_edges_description = hexlib.create_triangle_edges(triangle, num_pillars, pillars_thickness)
hexlib.add_new_edges(new_edges_description, vertices, edges)

# Now, transform to square every vertex we have:
vertices = np.array(vertices)
transformation_matrix = np.matrix('1.0 -0.577350269189626; 0.0 1.154700538379251')
resulting_vertices = transformation_matrix * vertices.transpose()
vertices = np.asarray(resulting_vertices.transpose())

# Finally, transpose to origin and reflect through the diagonal 'y = -x'
vertices -= 1.0

reflected_vertices = vertices[:, [1, 0]] * (-1)
reflected_edges = np.array(edges) + len(vertices)

# saving vertices on hypotenuses
hypotenuse_nodes = np.array(hypotenuse_nodes)
resulting_hypotenuse_nodes = transformation_matrix * hypotenuse_nodes.transpose()
hypotenuse_nodes = np.asarray(resulting_hypotenuse_nodes.transpose())
hypotenuse_nodes -= 1.0
reflected_hypotenuse_nodes = hypotenuse_nodes[:, [1, 0]] * (-1)
hypotenuse_nodes = hypotenuse_nodes.tolist() + reflected_hypotenuse_nodes.tolist()

# saving vertices on pillars
pillar_nodes = np.array(pillar_nodes)
resulting_pillar_nodes = transformation_matrix * pillar_nodes.transpose()
pillar_nodes = np.asarray(resulting_pillar_nodes.transpose())
pillar_nodes -= 1.0
reflected_pillar_nodes = pillar_nodes[:, [1, 0]] * (-1)
pillar_nodes = pillar_nodes.tolist() + reflected_pillar_nodes.tolist()

# transforming pillar example
pillar_example = np.array(pillar_example)
resulting_pillar_example = transformation_matrix * pillar_example.transpose()
pillar_example = np.asarray(resulting_pillar_example.transpose())
delta_y = abs(pillar_example[1][1] - pillar_example[0][1])
norm_example = np.linalg.norm(pillar_example[1] - pillar_example[0])
thickness_correction_factor = delta_y / norm_example
#print("Thickness correction factor: " + str(thickness_correction_factor))

new_edges_description = []
for edge in edges:
    v1 = edge[0]
    v2 = edge[1]

    reflected_vertex1 = reflected_vertices[v1]
    reflected_vertex2 = reflected_vertices[v2]

    new_edges_description.append([reflected_vertex1, reflected_vertex2])

vertices = vertices.tolist()
hexlib.add_new_edges(new_edges_description, vertices, edges)

# triangle is also transformed and reflected
triangle = np.array(triangle)
resulting_triangle = transformation_matrix * triangle.transpose()
transformed_triangle = np.asarray(resulting_triangle.transpose())

transformed_triangle -= 1.0
reflected_triangle = transformed_triangle[:, [1, 0]] * (-1)

triangle_vertices = list(transformed_triangle) + list(reflected_triangle)

incenter_triangle_pairs = hexlib.add_polygons_incenters([transformed_triangle, reflected_triangle], vertices, edges)

# print wire output
hexlib.create_wire(vertices, edges, out_wire)


print("Inflating ...")

# Computing void thickness and necessary resolution
if num_pillars == 1:
    thickness_void = 0.0
    min_resolution = 2 / (2*thickness)
    chosen_resolution = math.pow(2, math.ceil(math.log(min_resolution) / math.log(2)))
else:
    thickness_void = (triangle_side - 2 * num_pillars * thickness) / (num_pillars - 1)
    min_resolution = max(2 / thickness_void, 2 / (2*thickness))
    chosen_resolution = math.pow(2, math.ceil(math.log(min_resolution) / math.log(2)))

if chosen_resolution > 2048:
        print("Resolution of " + str(chosen_resolution) + " is too big")
        print("Skipping experiment!")
        exit()

if chosen_resolution < 64:
    chosen_resolution = 64

print("Thickness void: " + str(thickness_void))
print("Minimum resolution: " + str(min_resolution))
print("Chosen resolution: " + str(chosen_resolution))

inflate_hexagonal_box_smarter(out_wire, str(0.00001), str(0.0000), out_mesh,
                              incenter_triangle_pairs + [[hypotenuse_nodes, float(thickness) * math.sqrt(2)],
                                                         [pillar_nodes, float(thickness)*thickness_correction_factor]], chosen_resolution)
