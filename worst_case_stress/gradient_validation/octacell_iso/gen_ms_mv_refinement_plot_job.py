#!/usr/bin/env python
import os, json, data
from collections import defaultdict
WCDir='/home/fjp234/microstructures/worst_case_stress'
expDir='{WCDir}/gradient_validation/octacell_iso'.format(WCDir=WCDir)
resultRoot='/scratch/fjp234/gradient_validation/wcs/octacell_iso_ms_mvol'
params = [0,2,4,5,9]

# Compute JFull/WCS/JS/gradient ranges for each P, param
# and output gnuplot ylim commands in rangeDir
rangeDir = resultRoot + '/ranges'
if not os.path.exists(rangeDir): os.makedirs(rangeDir)
rangeFilePath = lambda P, param, quantity: '{rangeDir}/P{P}_param{param}.{quantity}.range'.format(rangeDir=rangeDir, P=P, param=param, quantity=q)
quantities = ['JFull', 'WCS', 'JS']

for P in range(1, 7):
    for param in params:
        minVal = defaultdict(lambda : float('inf'))
        maxVal = defaultdict(lambda : float('-inf'))
        for degree in [1, 2]:
            for max_vol in [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]:
                for ms_grid_size in [128,256,512,1024,2048,4096]:
                    path='{root}/d{deg}/ms_{ms}/vol_{vol}/P{P}_param{param}.txt'.format(root=resultRoot, deg=degree, vol=max_vol, ms=ms_grid_size, P=P, param=param)
                    d = data.extractDataTable(path)
                    if (len(d['paramValues']) == 0): continue
                    for q in quantities:
                        minVal[q] = min(minVal[q], d[q].min())
                        minVal['D' + q] = min(minVal['D' + q], d['D' + q].min(), d['centeredDiff' + q].min(), d['forwardDiff' + q].min())
                        maxVal[q] = max(maxVal[q], d[q].max())
                        maxVal['D' + q] = max(maxVal['D' + q], d['D' + q].max(), d['centeredDiff' + q].max(), d['forwardDiff' + q].max())
        for q in sum([[qq, 'D' + qq] for qq in quantities], []):
            rangeFile = file(rangeFilePath(P, param, q), 'w')
            rangeFile.write('set yrange [{}:{}];\n'.format(minVal[q], maxVal[q]))

cmds = []
for degree in [1, 2]:
    for max_vol in [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]:
        for ms_grid_size in [128,256,512,1024,2048,4096]:
            resultDir='{root}/d{deg}/ms_{ms}/vol_{vol}'.format(root=resultRoot, deg=degree, vol=max_vol, ms=ms_grid_size)
            for P in range(1, 7):
                for param in params:
                    cmds.append(
                    {'cmd': "gnuplot -e \"resultDir='{resultDir}'; P={P}; paramNum={param}; rangeDir='{rangeDir}'\" plot.gpi".format(resultDir=resultDir, P=P, param=param, rangeDir=rangeDir),
                     'cwd': expDir});
print json.dumps(cmds, indent=2)
