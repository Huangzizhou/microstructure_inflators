#!/usr/bin/env python
import numpy as np;
from subprocess import check_call
import json

config = {
        "lp": 2,
        "radius": 0.4,
        "width": 10,
        "height": 10,
        "grid_size": [10, 10],
        "output": "tmp.obj"
        }

radius_range = np.arange(0.05, 0.5, 0.05);
p_range = np.arange(2.0, 10.0, 0.5);

for r in radius_range:
    for p in p_range:
        basename = "l{:0.3}_r{:0.3}_10x10".format(p, r);
        config["lp"] = p;
        config["radius"] = r;
        config["output"] = basename + ".obj";
        json_file = basename + ".config";
        with open(json_file, 'w') as fout:
            json.dump(config, fout);
        command = "../generate_holy_rectangle.py {}".format(json_file);
        check_call(command.split());


