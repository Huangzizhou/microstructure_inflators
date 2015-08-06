import math, numpy as np;
from math import log;
from numpy.linalg import norm

def log2(a): return log(a, 2)
def spacing(ERange, nuRange, nSubdiv):
    """
    Compute a spacing for E, nu in the passed ranges.
    E is evenly spaced on a logscale, and nuRange is evenly spaced on a
    linear scale. The goal is to use approximately nSubdiv samples for each,
    but to have sample point hit the whole (and fractional) powers of two.
    """
    youngInvExponent = math.ceil(nSubdiv / log2(ERange[1] / ERange[0]))
    poissonExponent = math.floor(log2((nuRange[1] - nuRange[0]) / nSubdiv))
    return (2 ** (1.0 / youngInvExponent), 2 ** poissonExponent)

def logYoungIterates(ERange, spacing):
    a = int(math.floor(log(ERange[0], spacing)))
    b = int(math.ceil(log(ERange[1], spacing)))
    return [spacing ** p for p in range(a, b + 1)] # a..b inclusive
    
def poissonIterates(nuRange, spacing):
    a = int(math.floor(nuRange[0] / spacing))
    b = int(math.ceil(nuRange[1] / spacing))
    return [spacing * p for p in range(a, b + 1)] # a..b inclusive

class Grid:
    """
    An infinite grid evenly spaced in (logE, nu) space. Assumes all Young's
    moduli are strictly positive.
    Also, bounds are specified on the Poisson's ratio to keep grid points
    feasible.
    """
    def __init__(self, logESpacing, nuSpacing, nuMin, nuMax):
        self.logESpacing = logESpacing
        self.nuSpacing = nuSpacing
        self.nuMin = nuMin
        self.nuMax = nuMax

    def pointAtIndices(self, Eidx, nuIdx):
        return (self.logESpacing ** Eidx, np.clip(self.nuSpacing * nuIdx, self.nuMin, self.nuMax))

    # Get the distance to a particular gridpoint in units of grid spacing.
    # For example, the distance to the closest gridpoint is always <= norm(0.5, 0.5)
    # Distance is in number of gridpoints
    def distToGridpoint(self, E, nu, gridEIdx, gridNuIdx):
        return norm(np.array(self.getPointIndex(E, nu, roundingOp = lambda x: x)) - (gridEIdx, gridNuIdx))

    # Find index of closest point
    def getPointIndex(self, E, nu, roundingOp = round): return map(roundingOp, (log(E, self.logESpacing), nu / self.nuSpacing))

    # Find closest gridpoint to (E, nu) in (logE, nu) space.
    def getPoint(self, E, nu, roundingOp = round): return pointAtIndices(*getPointIndex(E, nu, roundingOp))

    # Find closest gridpoint to (E, nu) in (logE, nu) space above or below
    def lowerBound(self, E, nu): return self.getPoint(E, nu, math.floor)
    def upperBound(self, E, nu): return self.getPoint(E, nu, math.ceil)
