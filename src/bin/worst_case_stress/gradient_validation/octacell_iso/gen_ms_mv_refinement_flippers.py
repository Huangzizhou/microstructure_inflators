#!/usr/bin/env python
import json, os
flippers = []
root='/scratch/fjp234/gradient_validation/wcs/octacell_iso_ms_mvol'
os.chdir(root)
for degree in [1, 2]:
    for max_vol in [1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6]:
        for ms_grid_size in [128,256,512,1024,2048,4096]:
            resultDir='d{deg}/ms_{ms}/vol_{vol}'.format(deg=degree, vol=max_vol, ms=ms_grid_size)
            frames = []
            for param in [0,2,4,5,9]:
                views = ['WCS', 'JS']
                frames = [{'image': ['P%i_param%i.%s.png' % (P, param, view) for view in views]} for P in range(1, 7)]
                flipper='{}/param{}.js'.format(resultDir, param)
                f = open(flipper, 'w')
                f.write("title = 'deg {deg} ms grid {ms} max_vol {mv} param {param} Gradient Validation';\n".format(deg=degree, ms=ms_grid_size, mv=max_vol, param=param))
                f.write("framesLabel = 'Global objective P:&nbsp;';\n");
                f.write("statistics = [];\n")
                f.write("views = [{}];\n".format(', '.join(["'%s'" % view for view in views])))
                f.write("frames = ")
                f.write(json.dumps(frames, indent=4))
                f.write(";\n")

                flippers.append(['param{param}:degree {deg}:boundary resolution {ms}:max_vol {mv:0.0e}'.format(param=param, deg=degree, ms=ms_grid_size, mv=max_vol), flipper])
f = open('directory.js', 'w')
f.write("flippers = %s;\n" % json.dumps(flippers, indent=4))
