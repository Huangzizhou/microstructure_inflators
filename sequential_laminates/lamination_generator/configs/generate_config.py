#!/usr/bin/env python
import argparse
import numpy as np
from mako.template import Template
import sys
import os.path

def parse_args():
    parser = argparse.ArgumentParser(description="Lamination config file generator");
    parser.add_argument("-r", "--rank", help="lamination rank", default=2,
            type=int);
    parser.add_argument("--ratio", help="material ratio",
            action="append", type=float);
    parser.add_argument("--angle", help="lamination angle in degrees",
            action="append", type=float);
    parser.add_argument("--scale", help="lamination refinement scale",
            action="append", type=float);
    args = parser.parse_args();
    return args;

def compute_total_material_ratio(ratios):
    ratio = ratios[0];
    for theta in ratios[1:]:
        ratio += (1.0 - ratio) * theta;
    return ratio;

def main():
    args = parse_args();
    assert(len(args.ratio) == args.rank);
    assert(len(args.angle) == args.rank);
    assert(len(args.scale) == args.rank);

    template_file = os.path.join(sys.path[0], "config.mako");
    config = Template(filename=template_file);

    total_material_ratio = compute_total_material_ratio(args.ratio);
    for num_triangle_per_unit in np.arange(1e5, 2e6, 1e5):
        num_triangle_per_unit = int(round(num_triangle_per_unit));
        ave_triangle_area = total_material_ratio / num_triangle_per_unit;
        contents = config.render(
                rank = args.rank,
                ratios = args.ratio,
                angles = args.angle,
                scales = args.scale,
                triangle_density = num_triangle_per_unit,
                max_area = ave_triangle_area
                );
        config_filename = "rank_2_laminate_{}.config".format(num_triangle_per_unit);
        with open(config_filename, 'w') as fout:
            fout.write(contents);

if __name__ == "__main__":
    main();
