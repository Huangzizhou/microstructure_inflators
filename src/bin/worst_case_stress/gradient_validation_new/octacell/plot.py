import subprocess, os, sys, numpy as np
import plot_gvalid

innerSmoothings = [0.009638, 0.01, 0.02]
internal_res = ["1e-4", "5e-5", "1e-5"]
bdry_res = [1024, 2048, 4096, 8192, 16384]

ismooth, = map(float, sys.argv[1:])
if ismooth not in innerSmoothings:
    raise Exception("Invalid innerSmoothing parameter")

MICRO_DIR=os.environ['MICRO_DIR']
WCDir=MICRO_DIR + '/worst_case_stress'
resultDir=os.environ['SCRATCH'] + '/gradient_validation/wcs/octacell_iso/ismooth{}'.format(ismooth)
if not os.path.exists(resultDir): raise Exception("Result directory not found: %s" % resultDir);

PValues = [1, 1.1, 2, 3, 4, 5, 6, 7, 8]
params = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13]
for P in PValues:
    for param in params:
        allWCS, allDWCS = np.zeros(0), np.zeros(0)
        for ir in internal_res:
            for br in bdry_res:
                pval, WCS, fdWCS, sdWCS = plot_gvalid.getData(resultDir + '/ires{}/bres{}/P{}_param{}.txt'.format(ir, br, P, param))
                allWCS = np.concatenate([allWCS, WCS])
                allDWCS = np.concatenate([allDWCS, fdWCS, sdWCS])
        # rangeWCS  = [0.9 * np.percentile(allWCS,  0.1), 1.1 * np.percentile(allWCS,  0.9)]
        # rangeDWCS = [0.9 * np.percentile(allDWCS, 0.1), 1.1 * np.percentile(allDWCS, 0.9)]
        rangeWCS  = [allWCS.min(), allWCS.max()];
        rangeDWCS = [allDWCS.min(), allDWCS.max()]

        print "P{} param{} WCS range: {}, DWCS range: {}".format(P, param, rangeWCS, rangeDWCS)
        continue

        for ir in internal_res:
            for br in bdry_res:
                pathPrefix = resultDir + '/ires{}/bres{}/P{}_param{}'.format(ir, br, P, param)
                plot_gvalid.makePlot(pathPrefix + '.txt', pathPrefix + '.png',
                                     rangeWCS, rangeDWCS)
                print "plotted '{}'".format(pathPrefix + '.png')
