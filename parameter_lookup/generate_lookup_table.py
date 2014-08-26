#!/usr/bin/env python

import argparse
import csv
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

def gather_modifier_files(modifier_dir):
    modifier_name_pattern = ".*\.modifier";
    return gather_files(modifier_dir, modifier_name_pattern);

def gather_result_files(result_dir):
    result_name_pattern = ".*_param\.json";
    return gather_files(result_dir, result_name_pattern);

def pair_up(modifier_files, result_files):
    modifier_matcher = re.compile("(.*)\.modifier");
    result_matcher = re.compile("(.*)_param\.json");

    modifier_basename = lambda f: modifier_matcher.match(os.path.basename(f)).group(1)
    result_basename = lambda f : result_matcher.match(os.path.basename(f)).group(1)

    modifier_dict = {modifier_basename(f):f for f in modifier_files};
    result_dict = {result_basename(f):f for f in result_files};

    keys = set();
    keys.update(modifier_dict.keys());
    keys.update(result_dict.keys());

    pairs = []
    for key in keys:
        if key in modifier_dict and key in result_dict:
            pairs.append((modifier_dict[key], result_dict[key]));
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
    for param in materials:
        parameters.append(param.values);
        if dim == 3:
            youngs.append([param.young_x, param.young_x, param.young_z]);
            shears.append([param.shear_yz, param.shear_zx, param.shear_xy]);
            poissons.append([
                param.poisson_yz, param.poisson_zy,
                param.poisson_zx, param.poisson_xz,
                param.poisson_xy, param.poisson_yz ]);
        elif dim == 2:
            youngs.append([param.young_x, param.young_x]);
            shears.append([param.shear_xy]);
            poissons.append([param.poisson_xy, param.poisson_yx ]);
    parameters = np.array(parameters);
    youngs = np.array(youngs);
    shears = np.array(shears);
    poissons = np.array(poissons);

    generate_and_save(parameters, os.path.join(out_dir, "all.index"));
    generate_and_save(youngs, os.path.join(out_dir, "young.index"));
    generate_and_save(shears, os.path.join(out_dir, "shear.index"));
    generate_and_save(poissons, os.path.join(out_dir, "poisson.index"));

    save_dataset(parameters, os.path.join(out_dir, "all.npy"));
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

def parse_args():
    parser = argparse.ArgumentParser(
            description="generate material to pattern parameter lookup table");
    parser.add_argument("--dim", type=int, help="dimension", choices=[2, 3],
            default=3);
    parser.add_argument("--output", "-o", help="output directory")
    parser.add_argument("modifier_dir", help="directory containing .modifier files");
    parser.add_argument("result_dir", help="directory containing _param.json files");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();

    modifier_files = gather_modifier_files(args.modifier_dir);
    result_files = gather_result_files(args.result_dir);
    file_pairs = pair_up(modifier_files, result_files);

    patterns = [];
    materials = [];
    for modifier_file, result_file in file_pairs:
        pattern_param = PatternParameter(modifier_file);
        material_param = MaterialParameter(args.dim, result_file);
        patterns.append(pattern_param);
        materials.append(material_param);

    generate_and_save_index(args.dim, materials, args.output);
    save_all_data(patterns, materials, args.output);
    save_data(materials, args.output, "material.csv");
    save_data(patterns, args.output, "pattern.csv");

if __name__ == "__main__":
    main();
