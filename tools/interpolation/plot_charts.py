#!/usr/bin/env python3
import re
import math
import argparse
import random

import plotly.graph_objects as go

import numpy as np

def plot2dCharts(original_values, optimized_values, axis_title='', title='', output="tmp_fig.html"):
    frames = list(range(1,len(original_values) + 1))

    trace_original = go.Scatter(
        name="Wasserstein",
        x=frames,
        y=original_values,
    )

    trace_optimized = go.Scatter(
        name="Optimized",
        x=frames,
        y=optimized_values,
    )

    data = [trace_original, trace_optimized]

    layout = go.Layout(
        title=title,
        yaxis=dict(
            autorange=True,
            title=axis_title,
            titlefont=dict(family='Courier New', size=18, color='#7f7f7f')
        ),
        xaxis=dict(
            title="Frames",
            titlefont=dict(family='Courier New', size=18, color='#7f7f7f')
        )

    )

    fig = go.Figure(data=data, layout=layout)
    fig.write_html(output, auto_open=True)



parser = argparse.ArgumentParser(description='Plot info about interpolation of microstructures.')

parser.add_argument('--stress-table',   help='table with worst case stresses')
parser.add_argument('--material-table', help='table with material properties')

args = parser.parse_args()

if args.stress_table:
    tableFile = open(args.stress_table)
    original_stresses = []
    optimized_stresses = []
    for line in tableFile:
        fields = line.strip().split()

        original_stresses.append(float(fields[1]))
        optimized_stresses.append(float(fields[2]))

    plot2dCharts(original_stresses, optimized_stresses, "WCS", "WCS", "wcs.html")

if args.material_table:
    tableFile = open(args.material_table)
    original_poisson = []
    optimized_poisson = []
    original_youngs = []
    optimized_youngs = []
    original_shear = []
    optimized_shear = []
    original_volume = []
    optimized_volume = []
    for line in tableFile:
        fields = line.strip().split()

        original_poisson.append(float(fields[1]))
        optimized_poisson.append(float(fields[2]))

        original_youngs.append(float(fields[3]))
        optimized_youngs.append(float(fields[4]))

        original_shear.append(float(fields[5]))
        optimized_shear.append(float(fields[6]))

        original_volume.append(float(fields[7]))
        optimized_volume.append(float(fields[8]))

    plot2dCharts(original_poisson, optimized_poisson, "Poisson's ratio", "Poisson's Ratio", "poisson.html")
    plot2dCharts(original_youngs, optimized_youngs, "Young's modulus", "Young's Modulus", "youngs.html")
    plot2dCharts(original_shear, optimized_shear, "Shear modulus", "Shear Modulus", "shear.html")
    plot2dCharts(original_volume, optimized_volume, "Volume fraction", "Volume Fraction", "volume.html")

