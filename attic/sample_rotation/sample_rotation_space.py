#!/usr/bin/env python
import argparse
import os.path

import numpy as np

from rotate_mesh import rotate, load_mesh, save_mesh

def parse_args():
    parser = argparse.ArgumentParser(description="Generate rotated models");
    parser.add_argument("input", help="input model");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    basename, ext = os.path.splitext(args.input);

    mesh = load_mesh(args.input);
    for theta in np.arange(0.0, 89.9, 10.0):
        for phi in np.arange(0.0, 89.9, 10.0):
            rot_mesh = rotate(mesh, theta, phi);
            out_name = "{}_theta_{:.1f}_phi_{:.1f}{}".format(
                    basename, theta, phi, ext);
            save_mesh(rot_mesh, out_name);

if __name__ == "__main__":
    main();
