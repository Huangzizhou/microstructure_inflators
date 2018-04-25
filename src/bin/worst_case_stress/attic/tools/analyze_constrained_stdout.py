#!/usr/bin/env python
# Print the iterate with lowest value of "statName" among those with tensor
# fitting objective value JS below JSThreshold
import numpy as np
import sys, re

opt_stdout,statName,JSThreshold = sys.argv[1:]

JSThreshold = float(JSThreshold)

class Iterate:
    def __init__(self, statName):
        self.__dict__ = {k: None for k in ["params", "moduli", "anisotropy", "JS", "stat"]}
        self.statName = statName
    def isComplete(self):
        for k, v in self.__dict__.items():
            if (v is None): return False
        return True
    def __str__(self):
        return "\n".join(["{}:\t{}".format(k, v)
            for k, v in self.__dict__.items() if (k != "statName") and (k !=
                "stat")] + ["{}:\t{}".format(self.statName, self.stat)])

incompleteIt = Exception("Incomplete iterate parsed")

iterates = []
it = None
for line in file(opt_stdout):
    if (line.find('p:\t') == 0):
        if not ((it is None) or it.isComplete()):
            print it
            raise incompleteIt
        it = Iterate(statName)
        # it.params = np.array(map(float, line.strip().split()[1:]))
        it.params = "\t".join(line.strip().split()[1:])
        iterates.append(it)
    # if (line.find('moduli:\t') == 0): it.moduli = np.array(map(float, line.strip().split()[1:]))
    if (line.find('moduli:\t') == 0): it.moduli = ", ".join(["%0.5f" % float(m) for m in line.strip().split()[1:]])
    if (line.find('anisotropy:\t') == 0): it.anisotropy = float(line.strip().split()[1])
    if (line.find('JS:\t') == 0): it.JS = float(line.strip().split()[1])
    if (line.find(statName + ':\t') == 0): it.stat = float(line.strip().split()[1])

if  not it.isComplete():
    print it;
    raise incompleteIt

iterates = filter(lambda ii: ii.JS < JSThreshold, iterates)
iterates.sort(lambda a, b: cmp(a.stat, b.stat))
print iterates[0]
