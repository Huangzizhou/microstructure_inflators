#!/usr/bin/env python
import os
import sys
import numpy as np
import hexlib
import matplotlib.pyplot as plt
import interactiveplotlib as ipl
import re
from subprocess import call



story_size = 10

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
    print "usage: ./pillar2quasipillar-storyline.py <size of triangles> <num of pillars> <thickeness of pillars> " \
          "<output folder> <output table>"
    print "example: ./pillar2quasipillar-storyline.py 1.0 3 0.3 quasi-pillar-instances story-line.txt"
    sys.exit(-1)

triangle_side = float(sys.argv[1])
num_pillars = int(sys.argv[2])
pillars_thickness = float(sys.argv[3])
folder_path = sys.argv[4]
table_path = sys.argv[5]

# create folder if it does not exist yet
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# create initial pillar structure
name_pillars = folder_path + '/hexagon-pillars-n{}-s{}-t{}'.format(num_pillars, triangle_side, pillars_thickness)
wire_name = name_pillars + '.wire'
mesh_name = name_pillars + '.msh'

if not os.path.isfile(mesh_name):
    cwd = os.getcwd()
    cmd = [cwd + '/hexa-many-pillars-creator.py', str(triangle_side), str(num_pillars), str(pillars_thickness),
       wire_name, mesh_name]
    call(cmd)

# now, create descendants
thickness_values = list(np.linspace(0.005, pillars_thickness, story_size))
thickness_values += [pillars_thickness + 0.005, pillars_thickness + 0.01]
for index, connection_thickness in enumerate(thickness_values):
    if index == (story_size-1):
        continue

    name_quasi_pillars = folder_path + '/hexagon-quasi-pillars-n{}-s{}-t{}-c{}'.format(num_pillars, triangle_side,
                                                                                       pillars_thickness,
                                                                                       connection_thickness)

    wire_name = name_quasi_pillars + '.wire'
    mesh_name = name_quasi_pillars + '.msh'

    if os.path.isfile(mesh_name):
        continue

    cwd = os.getcwd()
    cmd = [cwd + '/hexa-quasi-pillars-creator.py', str(triangle_side), str(num_pillars), str(pillars_thickness),
           str(connection_thickness),
           wire_name, mesh_name]
    call(cmd)

# now, run homogenization on all instances created
cwd = os.getcwd()
cmd = [cwd + '/run-homogenization.py', str(folder_path), str(table_path)]
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
    m = re.search('-c(.+?).msh', fields[5])
    if m:
        note = round(float(m.group(1)),3)

    annotes.append(note)
    p.append(fields[0])

leftVertex = (-1.0, 0.0)
topVertex = (0.35, 200)
rightVertex = (1.0, 0.0)

fig, ax = plt.subplots()

# draw line
num_instances = len(x)
x = np.array(x)
y = np.array(y)
fit = np.polyfit(x[1:num_instances-1], y[1:num_instances-1], deg=1)
ax.plot(x, fit[0] * x + fit[1], color='red')

ax.grid(True)
polygon = plt.Polygon([leftVertex, topVertex, rightVertex], fill=None, color='r')
ax.add_patch(polygon)

ax.set_xlabel(r'$\nu$')
ax.set_ylabel('E')

col = ax.scatter(x, y, c='red', marker='+', picker=True)

af = ipl.AnnoteFinder(x, y, annotes, ax=ax, color='r', start_showing=True)
fig.canvas.mpl_connect('button_press_event', af)

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

# line to kappa and to \mu
x3 = (fit[0] * x[1:num_instances-1] + fit[1]) / (2.*(1-x[1:num_instances-1]))
y3 = (fit[0] * x[1:num_instances-1] + fit[1]) / (2.*(1+x[1:num_instances-1]))
#x3 = 63.86 - 0.735 / (1-x[1:num_instances-1])

#ax2.plot(x3, y2[1:num_instances-1], color='red')
ax2.plot(x3, y3, color='b')

# ax.axis('auto')
ax2.set_xlabel('$\kappa$')
ax2.set_ylabel(r'$\mu$')

col = ax2.scatter(x2, y2, c='blue', marker='o', picker=True)

af = ipl.AnnoteFinder(x2, y2, annotes, ax=ax2, color='b', start_showing=True)
fig2.canvas.mpl_connect('button_press_event', af)

plt.show()