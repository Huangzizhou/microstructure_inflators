#!/usr/bin/env python
import math
import sys

from cycler import cycler
import numpy as np
import matplotlib.pyplot as plt
import interactiveplotlib as ipl
from matplotlib.pyplot import cm

import hexlib

num_lines = 10

fig, ax = plt.subplots()

ax.grid(True)

color=cm.rainbow(np.linspace(0,1,num_lines))
for i,c in zip(range(num_lines),color):
    triangle = hexlib.theoretical_triangle(200.0, 0.35, 1.0 * (i+1)/num_lines)
    polygon = plt.Polygon(triangle, fill=None, color=c)
    ax.add_patch(polygon)

ax.set_xlabel(r'$\nu$')
ax.set_ylabel('E')

ax.set_xlim([-1.1, 1.1])
ax.set_ylim([-0.1, 200.0])

fig2, ax2 = plt.subplots()

# plt.axes()
ax2.grid(True)

color=cm.rainbow(np.linspace(0,1,num_lines))
for i,c in zip(range(num_lines),color):
    rectangle = hexlib.theoretical_rectangle(200.0, 0.35, 1.0 * (i+1)/num_lines)
    polygon = plt.Polygon(rectangle, fill=None, color=c)
    ax2.add_patch(polygon)

ax2.set_xlabel('$\kappa$')
ax2.set_ylabel(r'$\mu$')

ax2.set_xlim([-0.1, 160.0])
ax2.set_ylim([-0.1, 80.0])

plt.show()
