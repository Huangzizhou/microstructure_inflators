#!/usr/bin/env python
import sys
import matplotlib as mpl
import pandas as pd
import numpy as np
mpl.use('agg')
from matplotlib import pyplot as plt

def getData(gvalidFile):
    gv_data = pd.read_table(gvalidFile)
    pval = gv_data['param']
    deltaP = pval[1] - pval[0]
    fdWCS = np.gradient(gv_data['WCS'], deltaP, edge_order=2)
    sdWCS = gv_data['gradp WCS']
    WCS   = gv_data['WCS']
    return (pval, WCS, fdWCS, sdWCS)

def makePlot(gvalidFile, outPath, rangeWCS = None, rangeDWCS = None):
    pval, WCS, fdWCS, sdWCS = getData(gvalidFile)

    plt.style.use('ggplot')
    mpl.rcParams.update({'font.size': 11})

    plt.figure(figsize=(12, 9))
    plt.grid(True)

    ax = plt.subplot(2, 1, 1)
    ax.set_xlabel('Param value')
    ax.set_ylabel('dWCS')
    plt.plot(pval, fdWCS, 'r-',
             pval, sdWCS, 'ko--');
    plt.ylim(rangeDWCS)
    ax.margins(0.05)

    ax = plt.subplot(2, 1, 2)
    ax.set_xlabel('Param value')
    ax.set_ylabel('WCS')
    plt.plot(pval, WCS, 'r-')
    plt.ylim(rangeWCS)
    ax.margins(0.05)

    plt.tight_layout(pad=1, h_pad=1);
    plt.savefig(outPath, dpi=100)
    plt.close()

if __name__ == "__main__":
    gvalidFile,gvalidPlot = sys.argv[1:]
    makePlot(gvalidFile, gvalidPlot)
