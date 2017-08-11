#!/usr/bin/env pythonw
import itertools
import json
import os
import random
import sys

import hexlib

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import interactiveplotlib as ipl

if len(sys.argv) < 4:
    print "usage: ./plotLutTable.py <type> <material> <table>"
    print "example: ./plotLutTable.py type default-material results_table1.txt results_table2.txt ..."
    sys.exit(-1)

dim = 2

materialPath = sys.argv[2]
if os.path.isfile(materialPath):
    print "Using material properties in: " + materialPath
else:
    materialPath = '../../materials/B9Creator.material'
    print "Using default properties in: " + materialPath

materialJson = json.load(file(materialPath))
chart_type = sys.argv[1]
baseE = materialJson['young']
baseNu = materialJson['poisson']

ax = plt.gca()
ax.grid(True)

# parsing information in Lookup tables and adding them to the chart
i = 0
x = []
y = []
annotes = []
alpha = []
p = []
for i in range(3, len(sys.argv)):
    tablePath = sys.argv[i]
    tableFile = open(tablePath)
    for line in tableFile:
        fields = line.strip().split()
        x.append(float(fields[2]))
        y.append(float(fields[1]))
        alpha.append(float(fields[3]))
        annotes.append(fields[5])
        p.append(fields[0])

initialColors = ["r", "b", "g", "m", "y", "c", "k", "grey"]
colorNames = []
for name, hex in mcolors.cnames.iteritems():
    colorNames.append(name)
random.shuffle(colorNames)
colorNames = itertools.chain(initialColors, colorNames)
colors = itertools.cycle(colorNames)

allPatterns = set(p)
colorMap = {}
legend_patches = []
for pattern in allPatterns:
    currentColor = next(colors)

    colorMap[pattern] = currentColor
    
    legendString = "pattern: " + str(pattern)
    legend_patches.append(Patch(color=currentColor, label=legendString))

z = []
for index, pattern in enumerate(p):
    color = colorMap[pattern]
    if alpha[index] > 1.0:
        alpha[index] = 1.0 / alpha[index]
    color_alpha = mcolors.to_rgba(color, alpha[index])
    z.append(color_alpha)
    #z.append(colorMap[pattern])

#print allPatterns
plt.close('all')
if chart_type == "triangle":
    leftVertex = (-1.0, 0.0)
    topVertex = (baseNu, baseE)
    rightVertex = (1.0, 0.0)

    fig, ax = plt.subplots()

    #plt.axes()
    ax.grid(True)
    polygon = plt.Polygon([leftVertex, topVertex, rightVertex], fill=None, color='r')
    ax.add_patch(polygon)

    #ax.axis('auto')
    ax.set_xlabel(r'$\nu$')
    ax.set_ylabel('E')

    col = ax.scatter(x, y, c=z, marker='+', picker=True)

    af = ipl.AnnoteFinder(x, y, annotes, ax=ax)
    fig.canvas.mpl_connect('button_press_event', af)

elif chart_type == "rectangle":

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

    fig, ax = plt.subplots()

    # plt.axes()
    ax.grid(True)

    # draw theoretical rectangle
    theoretical_rectangle = hexlib.theoretical_rectangle(200.0, 0.35, 1.0)
    polygon = plt.Polygon(theoretical_rectangle, fill=None, color='b')
    ax.add_patch(polygon)

    # ax.axis('auto')
    ax.set_xlabel('$\kappa$')
    ax.set_ylabel(r'$\mu$')

    col = ax.scatter(x2, y2, c=z, marker='o', picker=True)

    af = ipl.AnnoteFinder(x2, y2, annotes, ax=ax)
    fig.canvas.mpl_connect('button_press_event', af)

else:
    fig, ax = plt.subplots()

    ax.grid(True)
    triangle = hexlib.theoretical_triangle(baseE, baseNu, 1.0)
    polygon = plt.Polygon(triangle, fill=None, color='r')
    ax.add_patch(polygon)

    ax.set_xlabel(r'$\nu$')
    ax.set_ylabel('E')

    col = ax.scatter(x, y, c='red', marker='+', picker=True)

    af = ipl.AnnoteFinder(x, y, annotes, ax=ax, color='r', start_showing=False)
    fig.canvas.mpl_connect('button_press_event', af)

    # for index, note in enumerate(annotes):
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
    rectangle = hexlib.theoretical_rectangle(baseE, baseNu, 1.0)
    polygon = plt.Polygon(rectangle, fill=None, color='b')
    ax2.add_patch(polygon)

    # ax.axis('auto')
    ax2.set_xlabel('$\kappa$')
    ax2.set_ylabel(r'$\mu$')

    col = ax2.scatter(x2, y2, c='blue', marker='o', picker=True)

    af = ipl.AnnoteFinder(x2, y2, annotes, ax=ax2, color='b', start_showing=False)
    fig2.canvas.mpl_connect('button_press_event', af)

#plt.legend(handles=legend_patches)
plt.show()
