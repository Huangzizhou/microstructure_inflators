#!/usr/bin/env python

import argparse
import csv
import json
import os.path

from Parameter import Parameter

def load_result(param_file):
    with open(param_file, 'r') as fin:
        result = json.load(fin);
        return result;

def merge(parameters, results):
    assert(len(parameters) == len(results));
    data = {
            "young_x": [],
            "young_y": [],
            "young_z": [],
            "shear_xy": [],
            "shear_yz": [],
            "shear_zx": [],
            "poisson_yz": [],
            "poisson_zy": [],
            "poisson_zx": [],
            "poisson_xz": [],
            "poisson_xy": [],
            "poisson_yx": [],
            };
    for param, result in zip(parameters, results):
        names = param.names;
        values = param.values;
        for i,name in enumerate(names):
            if name in data:
                data[name].append(values[i]);
            else:
                data[name] = [values[i]];
        data["young_x"].append(result["youngs_modulus"][0]);
        data["young_y"].append(result["youngs_modulus"][1]);
        data["young_z"].append(result["youngs_modulus"][2]);
        data["shear_yz"].append(result["shear_modulus"][0]);
        data["shear_zx"].append(result["shear_modulus"][1]);
        data["shear_xy"].append(result["shear_modulus"][2]);
        data["poisson_yz"].append(result["poisson_ratio"][0]);
        data["poisson_zy"].append(result["poisson_ratio"][1]);
        data["poisson_zx"].append(result["poisson_ratio"][2]);
        data["poisson_xz"].append(result["poisson_ratio"][3]);
        data["poisson_xy"].append(result["poisson_ratio"][4]);
        data["poisson_yx"].append(result["poisson_ratio"][5]);
    return data;

def save_csv(data, output_file):
    basename, ext = os.path.splitext(output_file);
    csv_file = basename + ".csv";
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(data.keys());

        num_entries = len(data[data.keys()[0]]);
        for i in range(num_entries):
            row = [data[key][i] for key in data.keys()];
            writer.writerow(["{:6f}".format(entry) for entry in row]);

def save_arff(data, output_file):
    basename, ext = os.path.splitext(output_file);
    arff_file = basename + ".arff";
    with open(arff_file, 'w') as fout:
        fout.write("@relation parameter_sweep\n");
        for key in data.keys():
            fout.write("@attribute {} numeric\n".format(key));

        fout.write("@data\n");
        num_entries = len(data[data.keys()[0]]);
        for i in range(num_entries):
            row = [data[key][i] for key in data.keys()];
            fout.write(",".join(["{}".format(val) for val in row]) + "\n");

def parse_args():
    parser = argparse.ArgumentParser(
            description="prepare parameter sweep result for anslysis");
    parser.add_argument("--output", "-o", help="output file");
    parser.add_argument("param_files", help="fit result json files", nargs="+");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    parameters = [];
    results = [];
    for param_file in args.param_files:
        dirname, basename = os.path.split(param_file);
        name, ext = os.path.splitext(basename);
        modifier_file = name.replace("_param", "") + ".modifier";
        modifier_file = os.path.join("configs", modifier_file);
        assert(os.path.exists(modifier_file));

        parameters.append(Parameter(modifier_file));
        results.append(load_result(param_file));

    data = merge(parameters, results);
    save_csv(data, args.output);
    save_arff(data, args.output);

if __name__ == "__main__":
    main();
