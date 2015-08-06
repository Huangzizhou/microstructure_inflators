import numpy as np, re
from glob import glob

def reportStats(values):
    np.set_printoptions(linewidth=160)
    print "mean:\t",   np.mean(values, axis=0)
    print "min:\t",    np.min(values, axis=0)
    print "max:\t",    np.max(values, axis=0)
    print "median:\t", np.median(values, axis=0)
class LUT:
    # patternData: (pattern, youngNu, anisotropy, params)
    def __init__(self, db=None, patternData=None):
        if ((db == None) == (patternData == None)):
            raise Exception("Must specify either db or patternData")

        self.youngPoissonAnisotropy = []
        variableWidthParams = []
        self.pattern = []
        if db:
            for line in file(db):
                fields = line.split();
                self.pattern.append(int(fields[0]))
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
        # TODO: remove duplicates
        self.youngPoissonAnisotropy = np.concatenate((self.youngPoissonAnisotropy, b.youngPoissonAnisotropy))
        self.numParams = np.concatenate((self.numParams, b.numParams))
        self.params = np.concatenate((self.params, b.params))
        self.pattern = np.concatenate((self.pattern, b.pattern))

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
        self.youngPoissonAnisotropy = np.delete(self.youngPoissonAnisotropy, anisoIndices, 0)
        self.params                 = np.delete(self.params, anisoIndices, 0)
        self.pattern                = np.delete(self.pattern, anisoIndices, 0)

    def filterPattern(self, pat):
        notPatEntries = np.where(self.pattern != pat)
        self.youngPoissonAnisotropy = np.delete(self.youngPoissonAnisotropy, notPatEntries, 0)
        self.params                 = np.delete(self.params, notPatEntries, 0)
        self.pattern                = np.delete(self.pattern, notPatEntries, 0)

    def write(self, path):
        print "writing lut of size ", self.size(), " to ", path
        formatEntries = lambda l: ["%0.16f" % m for m in l]
        f = open(path, 'w')
        for i in range(self.size()):
            f.write("\t".join([str(self.pattern[i])] +
                formatEntries(self.youngPoissonAnisotropy[i]) + [str(len(self.params[i]))]
                              + formatEntries(self.params[i][0:self.numParams[i]])) + "\n")
        f.close()

def extract(pat, directory, printableOnly = True):
    # TODO: try also extracting only the last point of each run
    # (this is the one we want if patterm optimization is truly working)
    moduli, aniso, params, printable = [], [], [], []
    for fname in glob(directory + '/stdout_*.txt'):
        for line in file(fname):
            m = re.search('^moduli:\s*(\S*)\s*(\S*)\s*(\S*)\s*(\S*)', line)
            if (m): moduli.append(map(float, [m.group(1), m.group(3)]))
            p = re.search('^p:\s*(.*)', line)
            if (p): params.append(map(float, p.group(1).split()))
            a = re.search('^anisotropy:\s*(.*)', line)
            if (a): aniso.append(float(a.group(1)))
            prmatch = re.search('^printable:\s*(\S*)', line)
            if (prmatch): printable.append(prmatch.group(1) == '1')
        lengths = map(len, [moduli, aniso, params, printable])
        numIterates = min(lengths)
        if (numIterates != max(lengths)):
            sys.stderr.write("WARNING: invalid iterate printouts in '%s'\n" % fname)
            moduli    = moduli[:numIterates]
            aniso     = aniso[:numIterates]
            params    = params[:numIterates]
            printable = printable[:numIterates]
            continue
    # The first printed iterate is just the initial parameters
    moduli = moduli[1:]
    aniso  =  aniso[1:]
    params = params[1:]
    if printableOnly:
        # printability filter
        moduli = [m for m,p in zip(moduli,printable) if p]
        aniso  = [m for m,p in zip(aniso, printable) if p]
        params = [m for m,p in zip(params,printable) if p]
    return LUT(patternData=(pat, moduli, aniso, params))
