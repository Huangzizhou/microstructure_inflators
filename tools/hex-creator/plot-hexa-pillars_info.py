#!/usr/bin/env pythonw
import json
import os
import sys
import re

import matplotlib.pyplot as plt
import interactiveplotlib as ipl
import numpy as np

import hexlib


def plotCharts(info, colormap='viridis', colormap_legend='info', prefix = "temp-fig"):

    fig, ax = plt.subplots()
    ax.grid(True)
    triangle = hexlib.theoretical_triangle(baseE, baseNu, volume_fraction)
    polygon = plt.Polygon(triangle, fill=None, color='r')
    ax.add_patch(polygon)
    ax.set_xlabel(r'$\nu$')
    ax.set_ylabel('E')

    max_info = max(info)
    min_info = min(info)
    ticks = list(np.linspace(min_info, max_info, 10))

    my_cmap = plt.cm.get_cmap(colormap)
    col = ax.scatter(x, y, c=info, marker='+', cmap=my_cmap, picker=True)
    af = ipl.AnnoteFinder(x, y, annotes, ax=ax, color='r', start_showing=False)
    fig.canvas.mpl_connect('button_press_event', af)

    cbar = fig.colorbar(col, ticks=ticks)
    cbar.ax.set_ylabel(colormap_legend)

    fig.savefig(prefix + '-poissonXyoungs.png')

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
    rectangle = hexlib.theoretical_rectangle(baseE, baseNu, volume_fraction)
    polygon = plt.Polygon(rectangle, fill=None, color='b')
    ax2.add_patch(polygon)
    ax2.set_xlabel('$\kappa$')
    ax2.set_ylabel(r'$\mu$')

    col2 = ax2.scatter(x2, y2, c=info, cmap=my_cmap, marker='o', picker=True)
    af = ipl.AnnoteFinder(x2, y2, annotes, ax=ax2, color='b', start_showing=False)
    fig2.canvas.mpl_connect('button_press_event', af)

    cbar2 = fig2.colorbar(col2, ticks=ticks)
    cbar2.ax.set_ylabel(colormap_legend)
    
    fig2.savefig(prefix + '-bulkXyoungs.png')


if len(sys.argv) < 3:
    print "usage: ./plot-hexa-pillars-info.py <material> <table>"
    print "example: ./plot-hexa-pillars-info.py default-material results_table1.txt results_table2.txt ..."
    sys.exit(-1)

dim = 2
volume_fraction = 1.0

materialPath = sys.argv[1]
if os.path.isfile(materialPath):
    print "Using material properties in: " + materialPath
else:
    materialPath = '../../materials/Russia.material'
    print "Using default properties in: " + materialPath

materialJson = json.load(file(materialPath))
baseE = materialJson['young']
baseNu = materialJson['poisson']

#ax = plt.gca()
#ax.grid(True)

# parsing information in Lookup tables and adding them to the chart
i = 0
x = []
y = []
annotes = []
alpha = []
p = []
p1 = []
p2 = []
p3 = []
p4 = []
volfrac = []
for i in range(2, len(sys.argv)):
    tablePath = sys.argv[i]
    tableFile = open(tablePath)
    for line in tableFile:
        fields = line.strip().split()
        x.append(float(fields[2]))
        y.append(float(fields[1]))
        alpha.append(float(fields[3]))
        annotes.append(fields[5])
        p.append(fields[0])

        # parse parameters
        m = re.search('volfrac-(.+?)_', fields[5])
        if m:
            volfrac.append(float(m.group(1)))
        else:
            volfrac.append(0.0)
        
        m = re.search('p1-(.+?)_', fields[5])
        if m:
            p1.append(float(m.group(1)))
        else:
            p1.append(0.0)

        m = re.search('p2-(.+?)_', fields[5])
        if m:
            p2.append(float(m.group(1)))
        else:
            p2.append(0)

        m = re.search('p3-(.+?)_', fields[5])
        if m:
            p3.append(float(m.group(1)))
        else:
            p3.append(0)

        m = re.search('p4-(.+?).msh', fields[5])
        if m:
            p4.append(float(m.group(1)))
        else:
            p4.append(0)

plotCharts(volfrac, 'winter', 'volume fraction', 'volfrac-')
plotCharts(p1, 'winter', 'p1: triangle side', 'p1-')
plotCharts(p2, 'cool', 'p2: number of pillars', 'p2-')
plotCharts(p3, 'viridis', 'p3: chirality', 'p3-')
plotCharts(p4, 'spring', 'p4: thickness', 'p4-')

plt.show()
