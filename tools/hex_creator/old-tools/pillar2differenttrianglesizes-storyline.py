#!/usr/bin/env python
import os
import sys
import math
import numpy as np
import hexlib
import matplotlib.pyplot as plt
import interactiveplotlib as ipl
import re
from subprocess import call

story_size = 20

if len(sys.argv) != 7:
    print "usage: ./pillar2differenttrianglesizes-storyline.py <size of initial triangle> <size of last triangle> <num of pillars> <initial thickeness> " \
          "<output folder> <output table>"
    print "example: ./pillar2differenttrianglesizes-storyline.py 0.9 1.3 5 0.3 instances story-line.txt"
    sys.exit(-1)

initial_triangle_side = float(sys.argv[1])
final_triangle_side = float(sys.argv[2])
num_pillars = int(sys.argv[3])
initial_pillar_thickness = float(sys.argv[4])
folder_path = sys.argv[5]
table_path = sys.argv[6]

# create folder if it does not exist yet
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# create initial pillar structure
name_pillars = folder_path + '/hexagon-many-pillars-n{}-s{}-t{}'.format(num_pillars, initial_triangle_side, initial_pillar_thickness)
wire_name = name_pillars + '.wire'
mesh_name = name_pillars + '.msh'

# computing area, to be able to recompute what should be thickness when changing triangle size
h = math.sqrt(3)/3 * (1 - initial_triangle_side/2)
area = 3* (initial_pillar_thickness * num_pillars * h) + math.sqrt(3)/4 * initial_triangle_side**2
total_triangle_area = math.sqrt(3)/4 * 2**2
volume_fraction = area / total_triangle_area

print "Volume fraction: " + str(volume_fraction)

if not os.path.isfile(mesh_name):
    cwd = os.getcwd()
    cmd = [cwd + '/hexa-many-pillars-creator.py', str(initial_triangle_side), str(num_pillars), str(initial_pillar_thickness),
           wire_name, mesh_name]
    call(cmd)

# now, create descendants
triangle_side_values = np.linspace(initial_triangle_side, final_triangle_side, story_size)
for index, triangle_side in enumerate(triangle_side_values):
    new_h = math.sqrt(3)/3 * (1 - triangle_side/2)
    new_thickness = (area - math.sqrt(3)/4 * triangle_side**2) / (3 * num_pillars * new_h)

    if new_thickness < 1e-4:

        print "Thickness is negative =( !!"
        break;

    name_pillars = folder_path + '/hexagon-many-pillars-n{}-s{}-t{}'.format(num_pillars, triangle_side,
                                                                            new_thickness)

    wire_name = name_pillars + '.wire'
    mesh_name = name_pillars + '.msh'

    if os.path.isfile(mesh_name):
        continue

    cwd = os.getcwd()
    cmd = [cwd + '/hexa-many-pillars-creator.py', str(triangle_side), str(num_pillars), str(new_thickness),
           wire_name, mesh_name]
    call(cmd)

# now, run homogenization on all instances created
cwd = os.getcwd()
cmd = [cwd + '/run-homogenization.py', str(folder_path), str(table_path), '../../materials/Russia.material']
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
    m = re.search('-s(.+?)-t', fields[5])
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
