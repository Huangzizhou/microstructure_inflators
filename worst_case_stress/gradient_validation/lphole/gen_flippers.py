#!/usr/bin/env python
import json
flippers = []
for max_vol in ['1e-3', '3e-4', '1e-4', '3e-5', '1e-5']:
    for p_center in [1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8]:
        frames = []
        for param in range(2):
            views = ['WCS', 'JS']
            frames = [{'image': ['d1/vol_%s/pcenter_%0.1f/P%i_param%i.%s.png' % (float(max_vol), p_center, P, param, view) for view in views]} for P in range(1, 7)]
            flipper='vol{}_pcenter{}_param{}.js'.format(max_vol, p_center, param)
            f = open(flipper, 'w')
            f.write("title = 'max_vol {mv} p_center {pc} param {param} Gradient Validation';\n".format(mv=max_vol, pc=p_center, param=param))
            f.write("statistics = [];\n")
            f.write("views = [{}];\n".format(', '.join(["'%s'" % view for view in views])))
            f.write("frames = ")
            f.write(json.dumps(frames, indent=4))
            f.write(";\n")

            flippers.append(['p_center {pc}:max_vol {mv}:param{param}'.format(pc=p_center, mv=max_vol, param=param), flipper])
print "flippers = " + json.dumps(flippers, indent=4) + ";"
