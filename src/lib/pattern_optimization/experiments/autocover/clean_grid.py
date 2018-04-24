import sys, os
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import numpy as np, json, math
from numpy.linalg import norm

import LookupTable
from autocover_tools import materialSpaceGrid

lutPath, configPath, outPath = sys.argv[1:]

lut = LookupTable.LUT(lutPath)
config = json.load(file(configPath))
grid = materialSpaceGrid(config)

# compute current covering grid:
currERange = (lut.minE(), lut.maxE());
currNuRange = (lut.minNu(), lut.maxNu());
minIndices = map(int, grid.getPointIndex(currERange[0], currNuRange[0], math.floor))
maxIndices = map(int, grid.getPointIndex(currERange[1], currNuRange[1], math.ceil))


# Compute indices of all patterns
tableUnroundedIndices = np.zeros((lut.size(), 2))
for i, yna in enumerate(lut.youngPoissonAnisotropy):
    tableUnroundedIndices[i] = grid.getPointIndex(yna[0], yna[1], roundingOp=lambda x: x)
    closestPtIndices = map(round, tableUnroundedIndices[i])

# Choose the single closest pattern to each gridpoint
gridPoints = set()

deletePattern = np.ones(lut.size(), dtype=bool)

for EIdx in range(minIndices[0], maxIndices[0] + 1):
    for nuIdx in range(minIndices[1], maxIndices[1] + 1):
        targetE, targetNu = grid.pointAtIndices(EIdx, nuIdx)
        pt = np.array([EIdx, nuIdx])
        lutPointDists = norm(tableUnroundedIndices - pt, axis=1)
        minIdx = np.argmin(lutPointDists)
        # order = np.argsort(lutPointDists)

        # Skip gridpoints too far away from existing points
        if (lutPointDists[minIdx] > 2):
            continue

        deletePattern[minIdx] = False

lut.removeIndices(np.where(deletePattern))

lut.write(outPath)
