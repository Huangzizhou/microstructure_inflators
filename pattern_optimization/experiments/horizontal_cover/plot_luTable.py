#!/usr/bin/env pythonw
import json, itertools, sys, paths
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from LookupTable import LUT, extract as extractLUT
import shutil
import glob

if (len(sys.argv) < 3):
    print "usage: ./plotLutTable.py config table"
    print "example: ./plotLutTable.py 2D_autocover_config.json results_table1.txt results_table2.txt ..."
    sys.exit(-1)

# Drawing the theoretical triangle
configPath = sys.argv[1]
config = json.load(file(configPath))

dim = config['dim']
materialPath = paths.material(config['material'])
materialJson = json.load(file(materialPath))
baseE = materialJson['young']
baseNu = materialJson['poisson']

if (dim == 2):
    leftVertex = (-1.0, 0.0)
    topVertex = (baseNu, baseE)
    rightVertex = (1.0, 0.0)
elif (dim == 3):
    leftVertex = (-1.0, 0.0)
    topVertex = (baseNu, baseE)
    rightVertex = (0.5, 0.0)

plt.axes()
polygon = plt.Polygon([leftVertex, topVertex, rightVertex], fill=None, color='r')
plt.gca().add_patch(polygon)

ax = plt.gca()
ax.grid(True)
plt.axis('auto')
plt.xlabel(r'$\nu$')
plt.ylabel('E')

# parsing information in Lookup tables and adding them to the chart
i = 0
x = []
y = []
p = []
for i in range(2, len(sys.argv)):
    tablePath = sys.argv[i]
    tableFile = open(tablePath)
    for line in tableFile:
        fields = line.strip().split()
        x.append(float(fields[2]))
        y.append(float(fields[1]))
        p.append(float(fields[0]))

initialColors = ["r", "b", "g", "m", "y", "c", "k", "grey"]
colorNames = []
for name, hex in mcolors.cnames.iteritems():
    colorNames.append(name)
random.shuffle(colorNames)
colorNames = itertools.chain(initialColors, colorNames)
colors = itertools.cycle(colorNames)

allPatterns = sorted(set(p))
colorMap = {}
legend_patches = []
for pattern in allPatterns:
    currentColor = next(colors)

    colorMap[pattern] = currentColor
    
    legendString = "pattern " + str(pattern) 
    legend_patches.append(mpatches.Patch(color=currentColor, label=legendString))

z = []
for pattern in p:
    z.append(colorMap[pattern])

print allPatterns
#colors = cm.rainbow(z)
#colors = cm.autoscale(z)
plt.scatter(x, y, c=z, marker='+')
    
plt.legend(handles=legend_patches)
plt.show()
