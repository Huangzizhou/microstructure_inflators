#!/usr/bin/env python3
import sys

import hexlib

if len(sys.argv) != 4:
    print("usage: ./hexa-to-general-params.py <size of triangles (p1)> <num of pillars (p2)> <relative thickness (of pillars) (p4)>")
    print("example: ./hexa-to-general-params.py 0.8 3 0.9")
    sys.exit(-1)

p1 = float(sys.argv[1])
p2 = int(sys.argv[2])
p4 = float(sys.argv[3])

vertices, edges, custom_pairs = hexlib.generate_positive_poisson_topology_and_thickness_info(p1, p2, p4)

original_parameters = hexlib.generate_vertices_parameters(vertices, 0.00001, 0.0, custom_pairs)

parameters_string = ' '.join(str(param) for param in original_parameters)
print(parameters_string)


