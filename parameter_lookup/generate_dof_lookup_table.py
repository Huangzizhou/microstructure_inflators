#!/usr/bin/env python

import argparse
import csv
import numpy as np
import re
import os
import os.path

import pyflann

from MaterialParameter import MaterialParameter
from DofParameter import DofParameter

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

def generate_pattern_index_map(wire_list_file):
    wire_file_index_map = {};
    with open(wire_list_file, 'r') as fin:
        count = 0;
        for wire_file in fin:
            basename = os.path.basename(wire_file);
            name, ext = os.path.splitext(basename);
            assert(name not in wire_file_index_map);
            wire_file_index_map[name] = count;
            count += 1;

    return wire_file_index_map;

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

def save_all_data(pattern_index_map, patterns, materials, out_dir):
    pattern_names = [];
    pattern_indices = [];
    pattern_num_dofs = [];
    for pattern in patterns:
        name = pattern.pattern_name;
        pattern_names.append(name);
        pattern_indices.append(pattern_index_map[name]);
        pattern_num_dofs.append(len(pattern.dofs));

    max_num_dofs = np.amax(pattern_num_dofs);
    num_entries = len(materials);

    csv_file = os.path.join(out_dir, "fit.csv");
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        names = ["pattern_name", "pattern_id", "num_dofs"];
        names += ["dof_{}".format(i) for i in range(max_num_dofs)];
        names += materials[0].names;
        writer.writerow(names);
        for i in range(num_entries):
            pattern_name = pattern_names[i];
            pattern_index = pattern_indices[i];
            num_dofs = pattern_num_dofs[i];
            dofs = np.zeros(max_num_dofs);
            dofs[:num_dofs] = patterns[i].dofs;

            values = [pattern_name, pattern_index, num_dofs];
            values += dofs.tolist();
            values += materials[i].values;
            for i,entry in enumerate(values):
                if isinstance(entry, float):
                    values[i] = "{:6f}".format(entry);
            writer.writerow(values);

def save_data(properties, out_dir, name):
    csv_file = os.path.join(out_dir, name);
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(properties[0].names);
        for param in properties:
            writer.writerow(["{:6f}".format(val) for val in param.values]);

def save_pattern_parameters(patterns, pattern_index_map, out_dir, csv_file):
    pattern_names = [];
    pattern_indices = [];
    pattern_num_dofs = [];
    for pattern in patterns:
        name = pattern.pattern_name;
        pattern_names.append(name);
        pattern_indices.append(pattern_index_map[name]);
        pattern_num_dofs.append(len(pattern.dofs));
    max_num_dofs = np.amax(pattern_num_dofs);

    csv_file = os.path.join(out_dir, csv_file);
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        header = ["pattern_id"];
        header += ["dof_{}".format(i) for i in range(max_num_dofs)];
        writer.writerow(header);
        for i,pattern in enumerate(patterns):
            dofs = np.zeros(max_num_dofs);
            dofs[:pattern_num_dofs[i]] = pattern.dofs;
            values = [pattern_indices[i]] + \
                    ["{:6f}".format(val) for val in dofs];
            writer.writerow(values);

def print_summary(properties):
    if len(properties) == 0: return;
    names = properties[0].names;
    num_fields = len(names);

    values = np.array([p.values for p in properties]);
    assert(values.shape[1] == num_fields);
    min_values = np.amin(values, axis=0);
    max_values = np.amax(values, axis=0);

    for name, min_val, max_val in zip(names, min_values, max_values):
        print("{:>10} {:>10.3f} {:>10.3f}".format(name, min_val, max_val));

def parse_args():
    parser = argparse.ArgumentParser(
            description="generate material to pattern parameter lookup table");
    parser.add_argument("--dim", type=int, help="dimension", choices=[2, 3], default=3);
    parser.add_argument("--output", "-o", help="output index directory");
    parser.add_argument("wire_list_file", help="wire list file");
    parser.add_argument("config_dir", help="directory containing .config files");
    parser.add_argument("result_dir", help="directory containing _param.json files");
    return parser.parse_args();

def main():
    args = parse_args();
    pattern_index_map = generate_pattern_index_map(args.wire_list_file);
    config_files = gather_config_files(args.config_dir);
    file_pairs = pair_up(config_files, args.result_dir);

    pattern_params = [];
    material_params = [];
    for config_file, result_file in file_pairs:
        pattern_param = DofParameter(config_file);
        material_param = MaterialParameter(args.dim, result_file);
        pattern_params.append(pattern_param);
        material_params.append(material_param);

    generate_and_save_index(args.dim, material_params, args.output);
    save_all_data(pattern_index_map, pattern_params, material_params, args.output);
    save_data(material_params, args.output, "material.csv");
    save_pattern_parameters(pattern_params, pattern_index_map, args.output, "pattern.csv");
    print_summary(material_params);

if __name__ == "__main__":
    main();

