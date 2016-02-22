#!/usr/bin/env python
import json, os
flippers = []
root='/scratch/fjp234/gradient_validation/wcs/lphole_ns_mvol'
os.chdir(root)
for degree in [1, 2]:
    for max_vol in [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]:
        for nsubdiv in [128,512,2048,8192,16384]:
            for p_center in [1.2, 1.3, 1.4]:
                resultDir='d{deg}/ns_{ns}/vol_{vol}/pcenter_{pc}'.format(deg=degree, vol=max_vol, ns=nsubdiv, pc=p_center)
                frames = []
                for param in range(1):
                    views = ['WCS', 'JS']
                    frames = [{'image': ['P%i_param%i.%s.png' % (P, param, view) for view in views]} for P in range(1, 7)]
                    flipper='{}/param{}.js'.format(resultDir, param)
                    f = open(flipper, 'w')
                    f.write("title = 'deg {deg} nsubdiv {ns} max_vol {mv} param {param} LpHole Gradient Validation';\n".format(deg=degree, ns=nsubdiv, mv=max_vol, param=param))
                    f.write("framesLabel = 'Global objective P:&nbsp;';\n");
                    f.write("statistics = [];\n")
                    f.write("views = [{}];\n".format(', '.join(["'%s'" % view for view in views])))
                    f.write("frames = ")
                    f.write(json.dumps(frames, indent=4))
                    f.write(";\n")

                    flippers.append(['param{param}:p_center {pc}:degree {deg}:bdry segments {ns}:max_vol {mv:0.0e}'.format(param=param, deg=degree, ns=nsubdiv, mv=max_vol, pc=p_center), flipper])
f = open('directory.js', 'w')
f.write("flippers = %s;\n" % json.dumps(flippers, indent=4))
