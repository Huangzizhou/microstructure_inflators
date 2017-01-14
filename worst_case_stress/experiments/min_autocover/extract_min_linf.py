import re, sys
from glob import glob
import numpy as np

fname = sys.argv[1]
relModulusThreshold, poissonThreshold = map(float, sys.argv[2:])

def relDist(val, truth):
    return abs(val - truth) / truth

# print "\t".join(["P", "LpStress", "LinfStress", "params"])

targetM = None
printable, params, youngDist, shearDist, poissonDist, LpStress, LinfStress = [], [], [], [], [], [], []
for line in file(fname):
    if (line[0:14] == 'Target moduli:'):
        m = line.strip().split('\t')[1:]
        if (len(m) != 9): raise Exception("Invalid target moduli line: " + line)
        targetM = map(float, m)
    
    m = re.search('^moduli:\s*\S*.*', line)
    if (m):
        m = m.group(0).strip().split('\t')[1:]
        if (len(m) != 9): raise Exception("Invalid iterate moduli line " + line)
        m = map(float, m)

        if (targetM == None): raise Exception("Target moduli not set")

        youngDist.append(max(relDist(m[0], targetM[0]),
                             relDist(m[1], targetM[1]),
                             relDist(m[2], targetM[2])))
        poissonDist.append(max(abs(m[3] - targetM[3]),
                               abs(m[4] - targetM[4]),
                               abs(m[5] - targetM[5])));
        shearDist.append(max(relDist(m[6], targetM[6]),
                             relDist(m[7], targetM[7]),
                             relDist(m[8], targetM[8])))

    p = re.search('^p:\s*(.*)', line)
    if (p): params.append(map(float, p.group(1).split()))

    prmatch = re.search('^printable:\s*(\S*)', line)
    if (prmatch): printable.append(prmatch.group(1) == '1')

    wcs = re.search('^Max Ptwise WCS:\s*(.*)', line)
    if (wcs): LinfStress.append(float(wcs.group(1)))

    lp = re.search('^WCS:\s*(.*)', line)
    if (lp): LpStress.append(float(lp.group(1)))

lengths = map(len, [printable, params, youngDist, shearDist, poissonDist, LpStress, LinfStress])
# print lengths
numIterates = min(lengths)
if (numIterates != max(lengths)):
    sys.stderr.write("WARNING: invalid iterate printouts in '%s'\n" % fname)
    printable    =   printable[:numIterates]
    params       =      params[:numIterates]
    youngDist    =   youngDist[:numIterates]
    shearDist    =   shearDist[:numIterates]
    poissonDist  = poissonDist[:numIterates]
    LpStress     =    LpStress[:numIterates]
    LinfStress   =  LinfStress[:numIterates]

is_valid = lambda i: (printable[i]
                    and (  youngDist[i] < relModulusThreshold)
                    and (  shearDist[i] < relModulusThreshold)
                    and (poissonDist[i] < poissonThreshold))

# Get the minimum valid stress level
validIndices  = filter(is_valid, range(numIterates))
minStressIdx  = validIndices[np.array(LpStress)[validIndices].argmin()]
minLinfIdx    = validIndices[np.array(LinfStress)[validIndices].argmin()] 
minLinfStress = LinfStress[minLinfIdx]

reduction = minLinfStress / LinfStress[0]

print "\t".join(map(str, [fname, reduction, LinfStress[0], minLinfStress] + params[minLinfIdx]))
