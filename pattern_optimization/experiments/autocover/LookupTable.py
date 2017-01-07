import numpy as np, re
from glob import glob

def reportStats(values):
    np.set_printoptions(linewidth=160)
    print "mean:\t",   np.mean(values, axis=0)
    print "min:\t",    np.min(values, axis=0)
    print "max:\t",    np.max(values, axis=0)
    print "median:\t", np.median(values, axis=0)
class LUT:
    # patternData: (pattern, young, nu, anisotropy, params)
    def __init__(self, db=None, patternData=None):
        if ((db == None) == (patternData == None)):
            raise Exception("Must specify either db or patternData")

        self.youngPoissonAnisotropy = []
        variableWidthParams = []
        self.pattern = []
        if db:
            for line in file(db):
                fields = line.split();
                patname = fields[0]
                try: patname = int(patname)
                except: pass
                self.pattern.append(patname)
                self.youngPoissonAnisotropy.append(map(float, fields[1:4]))
                variableWidthParams.append(map(float, fields[5:]))
        elif patternData:
            pattern, youngNu, anisotropy, variableWidthParams = patternData
            if (len(youngNu) != len(variableWidthParams) or len(youngNu) != len(anisotropy)):
                raise Exception("Mismatched LUT data size")
            self.youngPoissonAnisotropy = [yn + [a] for yn, a in zip(youngNu, anisotropy)]
            self.pattern = len(youngNu) * [pattern]

        self.youngPoissonAnisotropy = np.array(self.youngPoissonAnisotropy)
        self.numParams = np.array([len(p) for p in variableWidthParams])
        self.pattern   = np.array(self.pattern)

        # pad params arrays to the largest size
        self.params = np.zeros((self.size(), np.max(self.numParams)))
        for i, p in enumerate(variableWidthParams):
            self.params[i,0:self.numParams[i]] = p

    def size(self):  return self.youngPoissonAnisotropy.shape[0]

    def minE(self):  return np.min(self.youngPoissonAnisotropy[:, 0])
    def maxE(self):  return np.max(self.youngPoissonAnisotropy[:, 0])

    def minNu(self): return np.min(self.youngPoissonAnisotropy[:, 1])
    def maxNu(self): return np.max(self.youngPoissonAnisotropy[:, 1])

    def union(self, b):
        self.youngPoissonAnisotropy = np.concatenate((self.youngPoissonAnisotropy, b.youngPoissonAnisotropy))
        self.numParams = np.concatenate((self.numParams, b.numParams))
        self.params = np.concatenate((self.params, b.params))
        self.pattern = np.concatenate((self.pattern, b.pattern))
        self.removeDuplicates()

    def reportStats(self):
        pats = np.unique(self.pattern)
        print "patterns: ", " ".join(map(str, pats))

        if (len(pats) > 1):
            print "Overall Young, Poisson, Anisotropy stats:"
            reportStats(self.youngPoissonAnisotropy)

        for pat in pats:
            print "\nPATTERN %i" % pat
            mask = self.pattern == pat
            print "Young, Poisson, Anisotropy stats:"
            reportStats(self.youngPoissonAnisotropy[mask])
            print "Parameter stats:"
            nparam = np.unique(self.numParams[mask])
            if (len(nparam.shape) != 1): raise Exception("A pattern should have a fixed number of parameters")
            reportStats(self.params[mask, 0:nparam[0]])

    def filterAnisotropy(self, anisoMin, anisoMax):
        aniso = self.youngPoissonAnisotropy[:, 2]
        anisoIndices = np.where(np.logical_or(aniso < anisoMin, aniso > anisoMax))
        self.removeIndices(anisoIndices)

    def removeIndices(self, indices):
        self.youngPoissonAnisotropy = np.delete(self.youngPoissonAnisotropy, indices, 0)
        self.params                 = np.delete(self.params, indices, 0)
        self.pattern                = np.delete(self.pattern, indices, 0)

    def filterPattern(self, pat):
        notPatEntries = np.where(self.pattern != pat)
        self.youngPoissonAnisotropy = np.delete(self.youngPoissonAnisotropy, notPatEntries, 0)
        self.params                 = np.delete(self.params, notPatEntries, 0)
        self.pattern                = np.delete(self.pattern, notPatEntries, 0)

    def filterInvalidPoisson(self, poissonLimit = 0.5):
        bad = [i for i, yna in enumerate(self.youngPoissonAnisotropy) if yna[1] > poissonLimit]
        self.removeIndices(bad)

    # Side effect: sorts the lookup table by Young's modulus
    def removeDuplicates(self):
        if self.size() < 2: return

        # Sort by Young's modulus
        order = np.argsort(self.youngPoissonAnisotropy[:, 0])
        self.youngPoissonAnisotropy = self.youngPoissonAnisotropy[order, :]
        self.numParams = self.numParams[order]
        self.params = self.params[order, :]
        self.pattern = self.pattern[order]

        ypaDiff = np.diff(self.youngPoissonAnisotropy, axis=0)
        paramDiff = np.diff(self.params, axis=0)
        patternDiff = np.diff(self.pattern)
        fullDiff = np.max(np.abs(ypaDiff), axis=1) + np.max(np.abs(paramDiff), axis=1) + np.max(np.abs(patternDiff))
        mask = np.ones(self.size(), 'bool')
        mask[1:] = fullDiff != 0

        self.youngPoissonAnisotropy = self.youngPoissonAnisotropy[mask, :]
        self.numParams = self.numParams[mask]
        self.params = self.params[mask, :]
        self.pattern = self.pattern[mask]

    def write(self, path):
        print "writing lut of size ", self.size(), " to ", path
        formatEntries = lambda l: ["%0.16f" % m for m in l]
        f = open(path, 'w')
        for i in range(self.size()):
            f.write("\t".join([str(self.pattern[i])] +
                formatEntries(self.youngPoissonAnisotropy[i]) + [str(len(self.params[i]))]
                              + formatEntries(self.params[i][0:self.numParams[i]])) + "\n")
        f.close()

# init: don't ignore the first iterate--we're initializing an autocover
def extract(dim, pat, directory, printableOnly = True, init = False):
    # TODO: try also extracting only the last point of each run
    # (this is the one we want if patterm optimization is truly working)
    # Note: we mark patterns with significantly differing young's moduli on
    # different axes as infinitely anisotropic (the anisotropy ratio doesn't
    # account for this)
    moduli, aniso, params, printable = [], [], [], []
    allYoung = []
    for fname in glob(directory + '/stdout_*.txt'):
        for line in file(fname):
            # print(line)
            m = re.search('^moduli:\s*\S*.*', line)
            if (m):
                m = m.group(0).strip().split('\t')
                # 2D: m = ['moduli', Ex, Ey, nu_xy, mu_xy]
                # 3D: m = ['moduli', Ex, Ey, Ez, nu_yz, ...]
                if dim == 2: moduli.append(map(float, [m[1], m[3]]))
                if dim == 3: moduli.append(map(float, [m[1], m[4]]))
                if dim == 2: allYoung.append(map(float, m[1:3]))
                if dim == 3: allYoung.append(map(float, m[1:4]))
            p = re.search('^p:\s*(.*)', line)
            if (p): params.append(map(float, p.group(1).split()))
            a = re.search('^anisotropy:\s*(.*)', line)
            if (a): aniso.append(float(a.group(1)))
            prmatch = re.search('^printable:\s*(\S*)', line)
            if (prmatch): printable.append(prmatch.group(1) == '1')
        lengths = map(len, [moduli, aniso, params, printable])
        # print lengths
        numIterates = min(lengths)
        if (numIterates != max(lengths)):
            sys.stderr.write("WARNING: invalid iterate printouts in '%s'\n" % fname)
            moduli    = moduli[:numIterates]
            aniso     = aniso[:numIterates]
            params    = params[:numIterates]
            printable = printable[:numIterates]
            continue
    # The first printed iterate is just the initial parameters
    # Ignore it unless we're initializing an autocover (init == True)
    if not init:
        moduli   =   moduli[1:]
        aniso    =    aniso[1:]
        params   =   params[1:]
        allYoung = allYoung[1:]
    # mark as infinitely anisotropic the patterns with > 4% difference in
    # Young moduli
    for i, y in enumerate(allYoung):
        if abs((max(y) - min(y)) / min(y)) > 0.04:
            aniso[i] = float('inf')
    if printableOnly:
        # printability filter
        moduli = [m for m,p in zip(moduli,printable) if p]
        aniso  = [a for a,p in zip(aniso, printable) if p]
        params = [x for x,p in zip(params,printable) if p]
    l = LUT(patternData=(pat, moduli, aniso, params))
    # NOTE: ceres appears to print the same iterate multiple times (at least
    # they appear to be the same to output precision)
    l.removeDuplicates()
    return l
