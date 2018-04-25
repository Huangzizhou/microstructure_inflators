import json, os
MICRO_DIR=os.environ['MICRO_DIR']
WCDir=MICRO_DIR + '/worst_case_stress'

params = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13]
internal_res = ["1e-4", "5e-5", "1e-5"]
bdry_res = [1024, 2048, 4096, 8192, 16384]
PValues = [1, 1.1, 2, 3, 4, 5, 6, 7, 8]

flippers = []

for ires in internal_res:
    for bres in bdry_res:
        resultDir=os.environ['SCRATCH'] + '/gradient_validation/wcs/octacell_iso/ires{ir}/bres{br}'.format(ir=ires, br=bres)
        if not os.path.exists(resultDir): raise Exception("Result directory not found: %s" % resultDir);

        frames=[]
        views = ["p%i" % p for p in params]
        for P in PValues:
            frames.append({'image': ['P{}_param{}.png'.format(P, param) for param in params]})

        flipper = open(resultDir + '/frames.js', 'w')
        flipper.write("title = 'Gradient Validation, Octacell (Isosurface) Internal resolution %s, Boundary resolution %i';\n" % (ires, bres))
        flipper.write("statistics = [];\n")
        flipper.write("views = %s;\n" % json.dumps(views))
        flipper.write("framesLabel = 'Global objective P';\n")
        flipper.write("frames = %s;\n" % json.dumps(frames))
        flippers.append(['IR {}:BR {}'.format(ires, bres),
                         'ires{}/bres{}/frames.js'.format(ires, bres)])

print "flippers = %s;" % json.dumps(flippers, indent=2)
