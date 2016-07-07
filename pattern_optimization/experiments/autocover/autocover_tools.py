from spacing import *
from optimization import PatternOptimization
import paths, inflation, sys
import pattern_constraints
from LookupTable import LUT, extract as extractLUT
import os
import shutil, archive

# The grid is chosen so that a (~targetNSubdiv x ~NSubdivsubdiv) subset of
# lattice points fully cover the approxERange and approxNuRange ranges.
def materialSpaceGrid(config):
    approxERange =  config['targetERange']
    approxNuRange = config['targetNuRange']
    if (len(approxERange) != 2 or len(approxNuRange) != 2):
        raise Exception("Invalid approximate material property ranges")
    logESpace, nuSpacing = spacing(approxERange, approxNuRange, config['targetNSubdiv'])
    # Constrain nu to a range where the elasticity tensor is positive definite
    nuMin = -1.0
    nuMax =  0.5 if config['dim'] == 3 else 1.0
    return Grid(logESpace, nuSpacing, nuMin + 1e-3, nuMax - 1e-3)

# Create a PatternOptimization instance filled with jobs to expand the 
# current frontier (and fill in holes) of a given lookup table.
def coverageExpansionOptimizer(config, lut, grid, constraints):
    # within 1/10 of the distance between gridpoints counts as a hit
    hitThreshold = 0.1

    # compute current covering grid:
    currERange = (lut.minE(), lut.maxE());
    currNuRange = (lut.minNu(), lut.maxNu());
    minIndices = map(int, grid.getPointIndex(currERange[0], currNuRange[0], math.floor))
    maxIndices = map(int, grid.getPointIndex(currERange[1], currNuRange[1], math.ceil))

    # Compute adjacent gridpoints
    # loop from currIdx - 1 to currIdx + 1
    # (margin to expand range)
    gridPoints = set()
    for EIdx in range(minIndices[0] - 1, maxIndices[0] + 2):
        for nuIdx in range(minIndices[1] - 1, maxIndices[1] + 2):
            gridPoints.add((EIdx, nuIdx))
        
    # Discard gridpoints reached
    hitPoints = set()
    tableUnroundedIndices = np.zeros((lut.size(), 2))
    for i, yna in enumerate(lut.youngPoissonAnisotropy):
        tableUnroundedIndices[i] = grid.getPointIndex(yna[0], yna[1], roundingOp=lambda x: x)
        closestPtIndices = map(round, tableUnroundedIndices[i])
        if grid.distToGridpoint(yna[0], yna[1], *closestPtIndices) < hitThreshold:
            hitPoints.add(tuple(closestPtIndices))

    # hit points should be a proper subset of the grid we laid down
    if (not hitPoints < gridPoints): raise Exception("LUT hit outside grid")

    unreached = gridPoints.difference(hitPoints)

    opt = PatternOptimization(config)
    opt.setConstraints(constraints)
    for pt in unreached:
        targetE, targetNu = grid.pointAtIndices(*pt)
        lutPointDists = norm(tableUnroundedIndices - pt, axis=1)
        order = np.argsort(lutPointDists)
        # Skip gridpoints too far away from existing points
        # (so coverage grows outword more smoothly, and we don't waste time
        # trying to hit unreachable regions)
        if (lutPointDists[order[0]] > 2):
            continue

        # TODO: take N closest points that are within some minimal distance of each other.
        startLUTIndices = np.unique([order[0],
                                     order[min(math.ceil(0.05 * len(order)), len(order) - 1)],
                                     order[min(math.ceil(0.10 * len(order)), len(order) - 1)]])
        for i in startLUTIndices: opt.enqueueJob(targetE, targetNu, lut.params[i])
    return opt

def roundName(i): return "round_%04i" % i
def roundDirectory(i): return roundName(i)
def roundLUTPath(i): return roundName(i) + '.txt'

def analyzeRuns(dim, prevLUT, num, pat, printableOnly, isotropicOnly):
    if prevLUT == None:
        # the previous round's lookup table better exist (we need to union it)
        prevLUT = LUT(roundLUTPath(num - 1))
    rdir = roundDirectory(num)
    lut = extractLUT(dim, pat, rdir, printableOnly = printableOnly)

    # Archive the runs.
    archive = tarfile.open(rdir + '.tgz', 'w')
    archive.add(rdir)
    archive.close()
    shutil.rmtree(rdir)

    if isotropicOnly: lut.filterAnisotropy(0.95, 1.05)
    # TODO: determine convergence by diffing against previous round lookup table?
    # also, could remove duplicates by doing a sort | uniq
    lut.union(prevLUT)
    return lut

# Construct the optimizer for this autocover round (with jobs enqueued)
def autocoverRoundOptimizer(num, config):
    if (num < 1): raise Exception("Autocover rounds should be numbered 1, ...")
    prev = num - 1

    dim = config['dim']
    pat = config['pattern']
    printableOnly = config.get('printable', True)

    grid = materialSpaceGrid(config)
    constraints = pattern_constraints.lookup(pat, dim)
    if (not printableOnly): constraints = [] # no equality constraints for unprintable!

    # Read in previous round's lookup table, extracting it from finished runs
    # if it doesn't exist
    lut = None
    if (os.path.exists(roundLUTPath(prev))):
        lut = LUT(roundLUTPath(prev))
    elif (os.path.exists(roundName(prev))):
        lut = analyzeRuns(dim, None, prev, pat, printableOnly, config.get('isotropic', True))
        lut.write(roundLUTPath(prev))
    else: raise Exception("Previous round (%i) does not exist" % (num - 1))

    lut = LUT(roundLUTPath(num - 1))
    roundDir = roundName(num)
    return coverageExpansionOptimizer(config, lut, grid, constraints)
