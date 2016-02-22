#!/usr/bin/env python
# Parse shape velocity/gradient validation output and print table
#   paramValue JFull ShepDerivativeJFull centeredDiffJFull forwardDiffJFull wcObjective wcDirectDObj wcAdjointDObj wcCenteredDiffObj wcForwardDiffObj wcDirectKLTerm wcDirectAdvectTerm 
# assumes the parameter values are evenly spaced.
import re, numpy as np, math
from collections import defaultdict

def extractDataTable(dataPath):
    d = defaultdict(list);
    for line in file(dataPath):
        if (re.match("[0-9]+\t[0-9]", line)):
            fields = line.strip().split()
            d['paramValues'].append(fields[1])
            d['JFull'      ].append(fields[2])
            d['DJFull'     ].append(fields[3])
            d['WCS'        ].append(fields[4])
            d['DWCS'       ].append(fields[5])
            d['DWCSDirect' ].append(fields[6])
            d['JS'         ].append(fields[7])
            d['DJS'        ].append(fields[8])

    d['paramValues'] = np.array(map(float, d['paramValues']))
    d['JFull'      ] = np.array(map(float, d['JFull']))
    d['WCS'        ] = np.array(map(float, d['WCS']))
    d['JS'         ] = np.array(map(float, d['JS']))

    d['DJFull'] = np.array(map(float, d['DJFull']))
    d['DWCS'  ] = np.array(map(float, d['DWCS']))
    d['DJS'   ] = np.array(map(float, d['DJS']))

    if (len(d['paramValues']) == 0): return d

    d['centeredDiffJFull'] = np.gradient(d['JFull'], d['paramValues'][1] - d['paramValues'][0])
    d['forwardDiffJFull' ] = np.diff(d['JFull']) / np.diff(d['paramValues'])

    d['centeredDiffWCS']   = np.gradient(d['WCS'], d['paramValues'][1] - d['paramValues'][0])
    d['forwardDiffWCS' ]   = np.diff(d['WCS']) / np.diff(d['paramValues'])

    d['centeredDiffJS']    = np.gradient(d['JS'], d['paramValues'][1] - d['paramValues'][0])
    d['forwardDiffJS' ]    = np.diff(d['JS']) / np.diff(d['paramValues'])
    return d

if __name__ == "__main__":
    import sys
    dataPath = sys.argv[1]
    d = extractDataTable(dataPath)
    print "\t".join(['PARAM',
                     'JFull', 'SD_JFull', 'CD_JFull', 'FD_JFull',
                     'WCS',   'SD_WCS',   'CD_WCS',   'FD_WCS'  , 'SD_WCS_DIRECT',
                     'JS',    'SD_JS',    'CD_JS',    'FD_JS']);
    for i in range(len(d['paramValues'])):
        print "\t".join(map(str, [d['paramValues'][i],
            d['JFull'][i], d['DJFull'][i], d['centeredDiffJFull'][i], "N/A" if i == 0 else d['forwardDiffJFull'][i - 1],
            d['WCS'  ][i], d['DWCS'  ][i], d['centeredDiffWCS'  ][i], "N/A" if i == 0 else d['forwardDiffWCS'  ][i - 1], d['DWCSDirect'][i],
            d['JS'   ][i], d['DJS'   ][i], d['centeredDiffJS'   ][i], "N/A" if i == 0 else d['forwardDiffJS'   ][i - 1]]))
