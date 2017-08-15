#!/usr/bin/env python
import math
import sys

import numpy as np

import hexlib
from hexlib import inflate_hexagonal_box_smarter

if len(sys.argv) != 6:
    print "usage: ./hex-creator.py <number of pillars per triangle> <size of triangles> <thickeness (diameter) of pillars> " \
          "<output wire> <output mesh>"
    print "example: ./hex-creator.py 3 0.2 0.01 output.wire output.msh"
    sys.exit(-1)

num_pillars = int(sys.argv[1])
triangle_side = float(sys.argv[2])
pillars_thickness = float(sys.argv[3]) / 2 # rest of the code assumes it is the radius
out_wire = sys.argv[4]
out_mesh = sys.argv[5]
parallelogram_side = 2.0
l = parallelogram_side / 2.0
s = triangle_side

# graph constructed on left triangle
vertices = []
edges = []

# define important vertices of parallelogram
origin = np.array([0, 0])
triangle_centroid = np.array([l, l * math.sqrt(3) / 3.0])

print "Constructing " + out_wire + " ..."

# First side: create pillars of the left side (counter-clock orientation)
#     /.
#    /  .
#   /.....
parallelogram_start = np.array([(l + s / 2.0) / 2.0, (l + s / 2.0) * math.sqrt(3) / 2.0])
parallelogram_end = np.array([(l - s / 2.0) / 2.0, (l - s / 2.0) * math.sqrt(3) / 2.0])

radius = s * math.sqrt(3) / 6
f = triangle_centroid - [0, radius]
triangle_start = f + [0, math.sqrt(3) * s / 2.0]
triangle_end = f - [s / 2.0, 0]

unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

parallelogram_start_offset = parallelogram_start + pillars_thickness * unit_vector
parallelogram_end_offset = parallelogram_end - pillars_thickness * unit_vector
triangle_start_offset = triangle_start + pillars_thickness * unit_vector
triangle_end_offset = triangle_end - pillars_thickness * unit_vector

# add new vertices and edges to current sets
new_edges_description = hexlib.create_pillars([parallelogram_start_offset, parallelogram_end_offset], [triangle_start_offset, triangle_end_offset], num_pillars)
new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
hexlib.add_new_edges(new_edges_description, vertices, edges)

# Second side: create pillars of the bottom side (counter-clock orientation)
#     ..
#    .  .
#   ______
parallelogram_start = np.array([l - s / 2.0, 0])
parallelogram_end = np.array([l + s / 2.0, 0])

triangle_start = f - [s / 2.0, 0]
triangle_end = f + [s / 2.0, 0]

unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

parallelogram_start_offset = parallelogram_start + pillars_thickness * unit_vector
parallelogram_end_offset = parallelogram_end - pillars_thickness * unit_vector
triangle_start_offset = triangle_start + pillars_thickness * unit_vector
triangle_end_offset = triangle_end - pillars_thickness * unit_vector

# add new vertices and edges to current sets
new_edges_description = hexlib.create_pillars([parallelogram_start_offset, parallelogram_end_offset], [triangle_start_offset, triangle_end_offset],
                                              num_pillars)
new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
hexlib.add_new_edges(new_edges_description, vertices, edges)

# Third side: create pillars of the right side (counter-clock orientation)
#     .\
#    .  \
#   .....\
parallelogram_start = np.array([2 * l - (l - s / 2.0) / 2, (l - s / 2.0) * math.sqrt(3) / 2])
parallelogram_end = np.array([2 * l - (l + s / 2.0) / 2, (l + s / 2.0) * math.sqrt(3) / 2])

triangle_start = f + [s / 2.0, 0]
triangle_end = f + [0, math.sqrt(3) * s / 2.0]

unit_vector = (parallelogram_end - parallelogram_start) / np.linalg.norm(parallelogram_end - parallelogram_start)
unit_vector_triangle = (triangle_end - triangle_start) / np.linalg.norm(triangle_end - triangle_start)

parallelogram_start_offset = parallelogram_start + pillars_thickness * unit_vector
parallelogram_end_offset = parallelogram_end - pillars_thickness * unit_vector
triangle_start_offset = triangle_start + pillars_thickness * unit_vector
triangle_end_offset = triangle_end - pillars_thickness * unit_vector

# add new vertices and edges to current sets
new_edges_description = hexlib.create_pillars([parallelogram_start_offset, parallelogram_end_offset], [triangle_start_offset, triangle_end_offset],
                                              num_pillars)
new_edges_description += [[triangle_start, triangle_start_offset], [triangle_end, triangle_end_offset]]
hexlib.add_new_edges(new_edges_description, vertices, edges)

# save nodes on hypotenuse
hypotenuse_nodes = []
for edge_description in new_edges_description:
    hypotenuse_nodes.append(edge_description[0])
    hypotenuse_nodes.append(edge_description[1])

# Now, adding edges on the triangle
triangle = [f + [0, math.sqrt(3) * s / 2.0], f - [s / 2.0, 0], f + [s / 2.0, 0]]
new_edges_description = hexlib.create_triangle_edges(triangle, num_pillars, pillars_thickness)
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

print "Inflating ..."
inflate_hexagonal_box_smarter(out_wire, str(pillars_thickness), str(0.00001), out_mesh,
                              incenter_triangle_pairs + [[hypotenuse_nodes, float(pillars_thickness) * math.sqrt(2)], [triangle_vertices, 0.001]])