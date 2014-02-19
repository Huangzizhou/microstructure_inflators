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
    parser.add_argument("--with-top-bottom-plates", type=bool, default=True);
    parser.add_argument("--single-material", default=False, action="store_true");
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

    bool_to_str = lambda b: "true" if b else "false";

    template_file = os.path.join(sys.path[0], "config.mako");
    config = Template(filename=template_file);

    total_material_ratio = compute_total_material_ratio(args.ratio);
    for num_triangle_per_unit in np.arange(1e3, 2e4, 1e3):
        num_triangle_per_unit = int(round(num_triangle_per_unit));
        ave_triangle_area = total_material_ratio / num_triangle_per_unit;
        file_index = "{:06}".format(num_triangle_per_unit);
        contents = config.render(
                rank = args.rank,
                ratios = args.ratio,
                angles = args.angle,
                scales = args.scale,
                file_index = file_index,
                max_triangle_area = ave_triangle_area,
                with_top_bottom_plates = bool_to_str(args.with_top_bottom_plates),
                single_material = bool_to_str(args.single_material)
                );
        config_filename = "rank_2_laminates_{}.config".format(file_index);
        with open(config_filename, 'w') as fout:
            fout.write(contents);

if __name__ == "__main__":
    main();
