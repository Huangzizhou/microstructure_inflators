#!/usr/bin/env python

import argparse
import csv
import json
import os.path
import re

def get_keys():
    return [
            "p",
            "radius",
            "width",
            "height",
            "poisson_ratio_1",
            "poisson_ratio_2",
            "shear_modulus",
            "youngs_modulus_1",
            "youngs_modulus_2",
            "residual_error",
            "condition_number" ];

def extract_parameter_from_filename(filename):
    """ Expect filename to be generated by
    "l{}_r{}_{}x{}.param".format(p, radius, width, height);
    """
    basename = os.path.basename(filename);
    pattern = "l(\d+(\.\d+)?)_r(\d+(\.\d+)?)_(\d+)x(\d+).param";
    result = re.match(pattern, basename);

    p = float(result.group(1));
    radius = float(result.group(3));
    width = float(result.group(5));
    height = float(result.group(6));
    return p, radius, width, height;

def get_values(param_file):
    with open(param_file, 'r') as fin:
        param = json.load(fin);
        extra_param = extract_parameter_from_filename(param_file);

        return [
                extra_param[0],
                extra_param[1],
                extra_param[2],
                extra_param[3],
                param["poisson_ratio"][0],
                param["poisson_ratio"][1],
                param["shear_modulus"][0],
                param["youngs_modulus"][0],
                param["youngs_modulus"][1],
                param["residual_error"],
                param["condition_num"]
                ];

def param_to_csv(csv_file, param_files):
    with open(csv_file, 'w') as fout:
        csv_writer = csv.writer(fout);
        csv_writer.writerow(get_keys());

        for param_file in param_files:
            values = get_values(param_file);
            csv_writer.writerow(values);

def parse_args():
    parser = argparse.ArgumentParser(
            description="Combine all param file into csv file");
    parser.add_argument("-o", "--output", help="output csv file");
    parser.add_argument("param_files", nargs="*");
    return parser.parse_args();

def main():
    args = parse_args();
    param_to_csv(args.output, args.param_files);

if __name__ == "__main__":
    main();
