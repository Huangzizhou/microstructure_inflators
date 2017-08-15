#!/usr/bin/env python
import os
import sys
import math
import hexlib
import matplotlib.pyplot as plt
import interactiveplotlib as ipl
import re
from subprocess import call


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
    print "usage: ./singlepillar2manypillars-storyline.py <size of triangles> <max num of pillars> <thickeness of initial pillar> " \
          "<output folder> <output table>"
    print "example: ./singlepillar2manypillars-storyline.py 1.0 10 0.3 quasi-pillar-instances story-line.txt"
    sys.exit(-1)

triangle_side = float(sys.argv[1])
max_num_pillars = int(sys.argv[2])
initial_pillar_thickness = float(sys.argv[3])
folder_path = sys.argv[4]
table_path = sys.argv[5]

# create folder if it does not exist yet
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# compute area and volume fraction
h = math.sqrt(3)/3 * (1 - triangle_side/2)
area = 3* (initial_pillar_thickness * 1 * h) + math.sqrt(3)/4 * triangle_side**2
total_triangle_area = math.sqrt(3)/4 * 2**2
volume_fraction = area / total_triangle_area

print "Volume fraction: " + str(volume_fraction)

# create initial pillar structure
name_pillars = folder_path + '/hexagon-pillars-n{}-s{}-t{}'.format(1, triangle_side, initial_pillar_thickness)
wire_name = name_pillars + '.wire'
mesh_name = name_pillars + '.msh'

if not os.path.isfile(mesh_name):
    cwd = os.getcwd()
    cmd = [cwd + '/hexa-many-pillars-creator.py', str(triangle_side), str(1), str(initial_pillar_thickness),
           wire_name, mesh_name]
    call(cmd)

# now, create descendants
pillar_values = range(2, max_num_pillars + 1, 1)
for index, num_pillars in enumerate(pillar_values):
    pillars_thickness = initial_pillar_thickness / num_pillars

    name_pillars = folder_path + '/hexagon-many-pillars-n{}-s{}-t{}'.format(num_pillars, triangle_side,
                                                                            pillars_thickness)

    wire_name = name_pillars + '.wire'
    mesh_name = name_pillars + '.msh'

    if os.path.isfile(mesh_name):
        continue

    cwd = os.getcwd()
    cmd = [cwd + '/hexa-many-pillars-creator.py', str(triangle_side), str(num_pillars), str(pillars_thickness),
           wire_name, mesh_name]
    call(cmd)

# now, run homogenization on all instances created
cwd = os.getcwd()
cmd = [cwd + '/run-homogenization.py', str(folder_path), str(table_path), "../../materials/Russia.material"]
call(cmd)

# parse table file
i = 0
x = []
y = []
annotes = []
alpha = []
p = []
tableFile = open(table_path)
for line in tableFile:
    fields = line.strip().split()
    x.append(float(fields[2]))
    y.append(float(fields[1]))
    alpha.append(float(fields[3]))

    # find width of this table element
    note = fields[5]
    m = re.search('-n(.+?)-s', fields[5])
    if m:
        note = round(float(m.group(1)), 3)

    annotes.append(note)
    p.append(fields[0])

fig, ax = plt.subplots()

ax.grid(True)
triangle = hexlib.theoretical_triangle(1.0, 0.0, volume_fraction)
polygon = plt.Polygon(triangle, fill=None, color='r')
ax.add_patch(polygon)

ax.set_xlabel(r'$\nu$')
ax.set_ylabel('E')

col = ax.scatter(x, y, c='red', marker='+', picker=True)

af = ipl.AnnoteFinder(x, y, annotes, ax=ax, color='r', start_showing=True)
fig.canvas.mpl_connect('button_press_event', af)

#for index, note in enumerate(annotes):
#    ax.annotate(note, (x[index], y[index]))


# Now, create a new plot showing shear and bulk modulus
x2 = []
y2 = []
for index in range(0, len(x)):
    nu = x[index]
    E = y[index]

    K = E / (2 * (1 - nu))
    G = E / (2 * (1 + nu))

    x2.append(K)
    y2.append(G)

fig2, ax2 = plt.subplots()

# plt.axes()
ax2.grid(True)
rectangle = hexlib.theoretical_rectangle(1.0, 0.0, volume_fraction)
polygon = plt.Polygon(rectangle, fill=None, color='b')
ax2.add_patch(polygon)

# ax.axis('auto')
ax2.set_xlabel('$\kappa$')
ax2.set_ylabel(r'$\mu$')

col = ax2.scatter(x2, y2, c='blue', marker='o', picker=True)

af = ipl.AnnoteFinder(x2, y2, annotes, ax=ax2, color='b', start_showing=True)
fig2.canvas.mpl_connect('button_press_event', af)

plt.show()
