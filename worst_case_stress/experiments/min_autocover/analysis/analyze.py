import re, sys
from glob import glob
import numpy as np

relModulusThreshold, poissonThreshold = map(float, sys.argv[1:])

def relDist(val, truth):
    return abs(val - truth) / truth

print "\t".join(["path", "reduction", "initial", "optimized",
                 "targetE", "targetNu",
                 "minE", "minNu"])

# Note: only works in 3D for now.
for directory in glob('*_results'):
    for fname in glob(directory + '/stdout_*.txt'):
        targetM = None
        # distantYoung and distantPoisson: the x/y/z modulus with greatest
        # distance from the target.
        printable, params, youngDist, shearDist, poissonDist, distantYoung, distantPoisson, stress = [], [], [], [], [], [], [], []
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
                youngXYZDists = [relDist(m[0],  targetM[0]),
                                 relDist(m[1],  targetM[1]),
                                 relDist(m[2],  targetM[2])]
                poissonXYZDists = [  abs(m[3] - targetM[3]),
                                     abs(m[4] - targetM[4]),
                                     abs(m[5] - targetM[5])]
                shearXYZDists = [relDist(m[6],  targetM[6]),
                                 relDist(m[7],  targetM[7]),
                                 relDist(m[8],  targetM[8])]
                distantYoungIdx   = np.argmax(youngXYZDists)
                distantPoissonIdx = np.argmax(poissonXYZDists)
                distantShearIdx   = np.argmax(shearXYZDists)

                if (targetM == None): raise Exception("Target moduli not set")

                youngDist  .append(  youngXYZDists[distantYoungIdx])
                poissonDist.append(poissonXYZDists[distantPoissonIdx])
                shearDist  .append(  shearXYZDists[distantShearIdx])

                distantYoung  .append(m[distantYoungIdx])
                distantPoisson.append(m[3 + distantPoissonIdx])

            p = re.search('^p:\s*(.*)', line)
            if (p): params.append(map(float, p.group(1).split()))

            prmatch = re.search('^printable:\s*(\S*)', line)
            if (prmatch): printable.append(prmatch.group(1) == '1')

            wcs = re.search('^Max Ptwise WCS:\s*(.*)', line)
            if (wcs): stress.append(float(wcs.group(1)))

        lengths = map(len, [printable, params, youngDist, shearDist, poissonDist, stress])
        # print lengths
        numIterates = min(lengths)
        if (numIterates != max(lengths)):
            sys.stderr.write("WARNING: invalid iterate printouts in '%s'\n" % fname)
            printable   =   printable[:numIterates]
            params      =      params[:numIterates]
            youngDist   =   youngDist[:numIterates]
            shearDist   =   shearDist[:numIterates]
            poissonDist = poissonDist[:numIterates]
            stress      =      stress[:numIterates]
            continue

        try:
            is_valid = lambda i: (  printable[i]
                             and (  youngDist[i] < relModulusThreshold)
                             and (  shearDist[i] < relModulusThreshold)
                             and (poissonDist[i] < poissonThreshold))

            # Get the minimum valid stress level
            validIndices = filter(is_valid, range(numIterates))
            minStressIdx = validIndices[np.array(stress)[validIndices].argmin()]

            if (not is_valid(0)):
                # print printable[i], youngDist[i], shearDist[i], poissonDist[i]
                raise Exception("First iterate invalid")
            reduction = stress[minStressIdx] / stress[0]
            print "\t".join([fname, str(reduction), str(stress[0]), str(stress[minStressIdx]),
                             str(targetM[0]), str(targetM[3]),
                             str(distantYoung[minStressIdx]),
                             str(distantPoisson[minStressIdx])
                        ])
        except Exception as e:
            # print fname, "has no valid iterates (out of", numIterates,")"
            pass
