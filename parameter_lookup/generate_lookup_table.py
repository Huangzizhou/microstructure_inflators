#!/usr/bin/env python

import argparse
import csv
import json
import os
import os.path
import re
import numpy as np

import pyflann

from MaterialParameter import MaterialParameter
from PatternParameter import PatternParameter

def gather_files(dir_name, filename_pattern):
    assert(os.path.isdir(dir_name));
    matcher = re.compile(filename_pattern);

    files = [];
    for f in os.listdir(dir_name):
        r = matcher.match(f);
        if r is not None:
            files.append(os.path.join(dir_name, f));
    return files;

def gather_config_files(config_dir):
    config_name_pattern = ".*\.config";
    return gather_files(config_dir, config_name_pattern);

def gather_result_files(result_dir):
    result_name_pattern = ".*_param\.json";
    return gather_files(result_dir, result_name_pattern);

def pair_up(config_files, result_dir):
    pairs = [];
    for config_file in config_files:
        basename = os.path.splitext(os.path.basename(config_file))[0];
        result_file = os.path.join(result_dir, basename + "_param.json");
        if os.path.exists(result_file):
            pairs.append((config_file, result_file));
        else:
            print("Warning: {} is missing".format(result_file));
    return pairs;

def generate_and_save(data, output_name):
    table = pyflann.FLANN();
    table.build_index(data);
    table.save_index(output_name);

def save_dataset(data, output_name):
    np.save(output_name, data);

def generate_and_save_index(dim, materials, out_dir):
    parameters = [];
    youngs = [];
    shears = [];
    poissons = [];
    compliance_tensors = [];
    elasticity_tensors = [];
    for param in materials:
        parameters.append(param.values);
        compliance_tensors.append(param.compliance_tensor.ravel(order="C"));
        elasticity_tensors.append(param.elasticity_tensor.ravel(order="C"));
        if dim == 3:
            youngs.append([param.young_x, param.young_y, param.young_z]);
            shears.append([param.shear_yz, param.shear_zx, param.shear_xy]);
            poissons.append([
                param.poisson_yz, param.poisson_zy,
                param.poisson_zx, param.poisson_xz,
                param.poisson_xy, param.poisson_yz ]);
        elif dim == 2:
            youngs.append([param.young_x, param.young_y]);
            shears.append([param.shear_xy]);
            poissons.append([param.poisson_xy, param.poisson_yx ]);
    parameters = np.array(parameters);
    compliance_tensors = np.array(compliance_tensors);
    elasticity_tensors = np.array(elasticity_tensors);
    youngs = np.array(youngs);
    shears = np.array(shears);
    poissons = np.array(poissons);

    generate_and_save(parameters, os.path.join(out_dir, "all.index"));
    generate_and_save(compliance_tensors, os.path.join(out_dir, "compliance.index"));
    generate_and_save(elasticity_tensors, os.path.join(out_dir, "elasticity.index"));
    generate_and_save(youngs, os.path.join(out_dir, "young.index"));
    generate_and_save(shears, os.path.join(out_dir, "shear.index"));
    generate_and_save(poissons, os.path.join(out_dir, "poisson.index"));

    save_dataset(parameters, os.path.join(out_dir, "all.npy"));
    save_dataset(compliance_tensors, os.path.join(out_dir, "compliance.npy"));
    save_dataset(elasticity_tensors, os.path.join(out_dir, "elasticity.npy"));
    save_dataset(youngs, os.path.join(out_dir, "young.npy"));
    save_dataset(shears, os.path.join(out_dir, "shear.npy"));
    save_dataset(poissons, os.path.join(out_dir, "poisson.npy"));

def save_all_data(patterns, materials, out_dir):
    csv_file = os.path.join(out_dir, "fit.csv");
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        names = patterns[0].names + materials[0].names;
        writer.writerow(names);
        for p_param, m_param in zip(patterns, materials):
            values = p_param.values + m_param.values;
            writer.writerow(["{:6f}".format(val) for val in values]);

def save_data(properties, out_dir, name):
    csv_file = os.path.join(out_dir, name);
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(properties[0].names);
        for param in properties:
            writer.writerow(["{:6f}".format(val) for val in param.values]);

def save_modifier(config, out_dir):
    modifier_file = os.path.join(out_dir, "lookup.modifier");
    with open(modifier_file, 'w') as fout:
        json.dump(config, fout, indent=4);

def parse_args():
    parser = argparse.ArgumentParser(
            description="generate material to pattern parameter lookup table");
    parser.add_argument("--dim", type=int, help="dimension", choices=[2, 3],
            default=3);
    parser.add_argument("--output", "-o", help="output directory")
    parser.add_argument("config_dir", help="directory containing .config files");
    parser.add_argument("result_dir", help="directory containing _param.json files");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();

    config_files = gather_config_files(args.config_dir);
    file_pairs = pair_up(config_files, args.result_dir);

    if len(file_pairs) == 0:
        raise RuntimeError("No data extracted from inputs");

    patterns = [];
    materials = [];
    for config_file, result_file in file_pairs:
        pattern_param = PatternParameter(config_file);
        material_param = MaterialParameter(args.dim, result_file);
        patterns.append(pattern_param);
        materials.append(material_param);

    generate_and_save_index(args.dim, materials, args.output);
    save_all_data(patterns, materials, args.output);
    save_data(materials, args.output, "material.csv");
    save_data(patterns, args.output, "pattern.csv");
    save_modifier(patterns[0].modifier_config, args.output);

if __name__ == "__main__":
    main();
