#!/usr/bin/env python

import argparse
import csv
import numpy as np
import re
import os
import os.path

import pyflann

from IsotropicMaterial import IsotropicMaterial
from generate_dof_lookup_table import generate_and_save, save_dataset

def parse_sweep_summary_file(filename):
    pattern_ids = [];
    pattern_param = [];
    materials = [];
    with open(filename, 'r') as fin:
        line_count = 0;
        for row in fin:
            try:
                fields = row.split();
                assert(len(fields) > 5);
                pattern_id = int(fields[0]);
                young = float(fields[1]);
                poisson = float(fields[2]);
                anisotropy = float(fields[3]);
                num_dofs = int(fields[4]);
                dofs = [float(val) for val in fields[5:]];
                assert(len(dofs) >= num_dofs);
            except:
                print("Error parsing line {}".format(line_count));
                print(row);
                raise;

            if anisotropy > 1.2 or anisotropy < 5.0/6.0:
                continue;

            material = IsotropicMaterial(3, young, poisson);

            # Todo
            # Filter on printability

            pattern_ids.append(pattern_id);
            pattern_param.append(dofs);
            materials.append(material);
            line_count+=1;

    return pattern_ids, pattern_param, materials;

def generate_and_save_index(materials, out_dir):
    compliance_tensors = [];
    elasticity_tensors = [];
    for material in materials:
        compliance_tensors.append(material.compliance_tensor.ravel(order="C"));
        elasticity_tensors.append(material.elasticity_tensor.ravel(order="C"));

    compliance_tensors = np.array(compliance_tensors);
    elasticity_tensors = np.array(elasticity_tensors);

    generate_and_save(compliance_tensors, os.path.join(out_dir,
        "compliance.index"));
    generate_and_save(elasticity_tensors, os.path.join(out_dir,
        "elasticity.index"));
    save_dataset(compliance_tensors, os.path.join(out_dir, "compliance.npy"));
    save_dataset(elasticity_tensors, os.path.join(out_dir, "elasticity.npy"));

def get_pattern_names(pattern_ids):
    pattern_names = ["pattern{:04}".format(val) for val in pattern_ids];
    return pattern_names;

def get_pattern_dofs(pattern_params):
    return [len(dofs) for dofs in pattern_params];

def save_all_data(pattern_ids, pattern_params, materials, out_dir):
    pattern_names = get_pattern_names(pattern_ids);
    pattern_indices = pattern_ids;
    pattern_num_dofs = get_pattern_dofs(pattern_params);
    max_num_dofs= np.amax(pattern_num_dofs);
    num_entries = len(materials);

    csv_file = os.path.join(out_dir, "fit.csv");
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        names = ["pattern_name", "pattern_id", "num_dofs"];
        names += ["dof_{}".format(i) for i in range(max_num_dofs)];
        names += ["young", "poisson"];

        writer.writerow(names);
        for i in range(num_entries):
            pattern_name = pattern_names[i];
            pattern_index = pattern_indices[i];
            num_dofs = pattern_num_dofs[i];
            dofs = np.zeros(max_num_dofs);
            dofs[:num_dofs] = pattern_params[i];

            values = [pattern_name, pattern_index, num_dofs];
            values += dofs.tolist();
            values += [materials[i].young, materials[i].poisson];
            for i,entry in enumerate(values):
                if isinstance(entry, float):
                    values[i] = "{:6f}".format(entry);
            writer.writerow(values);

def save_material_parameters(materials, out_dir, out_name):
    csv_file = os.path.join(out_dir, out_name);
    header = ["young", "poisson"];
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(header);
        for material in materials:
            writer.writerow([material.young, material.poisson]);

def save_pattern_parameters(pattern_params, pattern_ids, out_dir, out_name):
    pattern_names = get_pattern_names(pattern_ids);
    pattern_indices = pattern_ids;
    pattern_num_dofs = get_pattern_dofs(pattern_params);
    max_num_dofs= np.amax(pattern_num_dofs);

    csv_file = os.path.join(out_dir, out_name);
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        header = ["pattern_id"];
        header += ["dof_{}".format(i) for i in range(max_num_dofs)];
        writer.writerow(header);

        for i, params in enumerate(pattern_params):
            dofs = np.zeros(max_num_dofs);
            dofs[:pattern_num_dofs[i]] = params;
            values = [pattern_indices[i]] + \
                    ["{:.6f}".format(val) for val in dofs];
            writer.writerow(values);

def parse_args():
    parser = argparse.ArgumentParser(
            description="generate material to pattern parameter lookup table");
    parser.add_argument("--output", "-o", help="output index directory");
    parser.add_argument("sweep_summary_file", help="sweep summary file");
    return parser.parse_args();

def main():
    args = parse_args();
    pattern_ids, pattern_params, materials =\
        parse_sweep_summary_file(args.sweep_summary_file);

    generate_and_save_index(materials, args.output);
    save_all_data(pattern_ids, pattern_params, materials, args.output);
    save_material_parameters(materials, args.output, "material.csv");
    save_pattern_parameters(pattern_params, pattern_ids, args.output, "pattern.csv");

if __name__ == "__main__":
    main();

