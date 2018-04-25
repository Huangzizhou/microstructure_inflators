#!/usr/bin/env python
import sys, re, numpy as np, math, pandas as pd
dataPath = sys.argv[1]
key = sys.argv[2]
if key not in {"pre-remesh", "post-remesh"}:
    raise Exception("Invalid data key");

data = {"wcs": [], "grad_wcs_norm": []}
for line in file(dataPath):
    m = re.match("WCS:\t(.+)", line)
    if (m): data["wcs"].append(float(m.group(1)))
    m = re.match(r"\|\|grad_p WCS\|\|:\t(.+)", line)
    if (m): data["grad_wcs_norm"].append(float(m.group(1)))

df = pd.DataFrame(data);
numIterates = len(df);
if (numIterates % 2 == 1): raise Exception("Error: odd number of iterates")
numIterates /= 2;
pre_remesh =  df[0:2 * numIterates:2]
post_remesh = df[1:2 * numIterates:2]

dataSet = pre_remesh if (key == "pre-remesh") else post_remesh

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
