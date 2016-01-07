#!/usr/bin/env python
import sys, re, numpy as np, math, pandas as pd
dataPath = sys.argv[1]

data = {"wcs": [], "grad_wcs_norm": [], "jvol": [], "grad_jvol_norm": [], "js": [], "grad_js_norm": []}
for line in file(dataPath):
    m = re.match("WCS:\t(.+)", line)
    if (m): data["wcs"].append(float(m.group(1)))
    m = re.match(r"\|\|grad_p WCS\|\|:\t(.+)", line)
    if (m): data["grad_wcs_norm"].append(float(m.group(1)))

    m = re.match("JS:\t(.+)", line)
    if (m): data["js"].append(float(m.group(1)))
    m = re.match(r"\|\|grad_p WCS\|\|:\t(.+)", line)
    if (m): data["grad_js_norm"].append(float(m.group(1)))

    m = re.match("JVol:\t(.+)", line)
    if (m): data["jvol"].append(float(m.group(1)))
    m = re.match(r"\|\|grad_p Jvol\|\|:\t(.+)", line)
    if (m): data["grad_jvol_norm"].append(float(m.group(1)))

dataSet = pd.DataFrame(data);
numIterates = len(dataSet);

print "i",
for c in dataSet.keys():
    print "\t",c,
print ""
for i in range(numIterates):
    print i,
    row = dataSet.iloc[i]
    for c in row.keys():
        print "\t",row[c],
    print ""
