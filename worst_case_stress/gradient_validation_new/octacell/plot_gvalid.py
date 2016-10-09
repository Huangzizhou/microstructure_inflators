#!/usr/bin/env python
import sys
import matplotlib as mpl
import pandas as pd
import numpy as np
mpl.use('agg')
from matplotlib import pyplot as plt

def getData(gvalidFile, stat):
    gv_data = pd.read_table(gvalidFile)
    pval = gv_data['param']
    deltaP = pval[1] - pval[0]
    fd= np.gradient(gv_data[stat], deltaP, edge_order=2)
    sd= gv_data['gradp ' + stat]
    val   = gv_data[stat]
    return (pval, val, fd, sd)

def makePlot(gvalidFile, outPath, stat, rangeVal = None, rangeD = None):
    pval, val, fd, sd= getData(gvalidFile, stat)

    plt.style.use('ggplot')
    mpl.rcParams.update({'font.size': 11})

    plt.figure(figsize=(12, 9))
    plt.grid(True)

    ax = plt.subplot(2, 1, 1)
    ax.set_xlabel('Param value')
    ax.set_ylabel('d' + stat)
    plt.plot(pval, fd, 'r-',
             pval, sd, 'ko--');
    plt.ylim(rangeD)
    ax.margins(0.05)

    ax = plt.subplot(2, 1, 2)
    ax.set_xlabel('Param value')
    ax.set_ylabel(stat)
    plt.plot(pval, val, 'r-')
    plt.ylim(rangeVal)
    ax.margins(0.05)

    plt.tight_layout(pad=1, h_pad=1);
    plt.savefig(outPath, dpi=100)
    plt.close()

if __name__ == "__main__":
    gvalidFile,gvalidPlot,stat = sys.argv[1:]
    makePlot(gvalidFile, gvalidPlot, stat)
