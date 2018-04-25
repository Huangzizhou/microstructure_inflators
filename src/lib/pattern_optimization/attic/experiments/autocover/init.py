from spacing import *
from homogenization import *
from optimization import PatternOptimization
import paths, inflation, sys
import pattern_constraints
from LookupTable import LUT

# E_min/E_max, nu_min/nu_max, and nsubdiv are for fitting a logspace grid spacing
if not (len(sys.argv) == 4 or len(sys.argv) == 5):
    print "usage: ./init.py dim pattern material ['init_params']"
    print "example: ./init.py 2 0 b9creator"
    sys.exit(-1)

maxRoundIters = 50

dim = int(sys.argv[1])
pat = int(sys.argv[2])
mat = paths.material(sys.argv[3])

def roundName(i):    return "round_%04i" % i
def roundLUTPath(i): return roundName(i) + '.txt'

constraints = pattern_constraints.lookup(pat, dim)

paramTypes = inflation.getParameterTypes(pat, dim, constraints)
initParams = None
if (len(sys.argv) == 5):
    initParams = map(float, sys.argv[4].strip().split())
    if (len(initParams) != len(paramTypes)):
        raise Exception("Invalid number of parameters")

if (initParams == None):
    # Default: middle of the thickness range, no offset
    def defaultParamValue(t):
        if (t == "thickness"): return 0.5 if dim == 3 else 0.1
        if (t == "offset"): return 0
        raise Exception("Unknown pattern type: " + t)
    initParams = [defaultParamValue(t) for t in paramTypes]

if (not inflation.isPrintable(pat, initParams, dim=dim, constraints=constraints)):
    raise Exception("Unprintable initial parameters")

young, poisson, aniso = homogenize_cubic(pat, initParams, mat, dim=dim, constraints=constraints)

print "Initing experiment at ", initParams
f = open(roundLUTPath(0), 'w')
f.write('\t'.join(map(str, [pat, young, poisson, aniso, len(initParams)] + initParams)))
f.write('\n')
f.close()
