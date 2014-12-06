#!/usr/bin/env python

import argparse
import csv
import os.path
import re
import numpy as np
from numpy.linalg import lstsq
from subprocess import check_call

def load_csv(csv_file):
    assert(os.path.exists(csv_file));
    with open(csv_file, 'r') as fin:
        reader = csv.reader(fin);
        header = reader.next();
        data = [];
        for row in reader:
            data.append(map(float, row));
        return header, np.array(data);

def select_columns(header, data, matchers):
    selected_columns = [];
    selected_headers = [];
    for i,name in enumerate(header):
        for f in matchers:
            if f(name) is not None:
                selected_headers.append(name);
                selected_columns.append(i);
                break;

    return selected_headers, data[:,selected_columns];

def extract_pattern_parameter(header, data):
    vertex_thickness_pattern = "vertex_orbit_\d+_thickness";
    edge_thickness_pattern = "edge_orbit_\d+_thickness";
    vertex_offset_pattern = "vertex_orbit_\d+_offset";
    matcher = [
            lambda s: re.search(vertex_thickness_pattern, s),
            lambda s: re.search(edge_thickness_pattern, s),
            lambda s: re.search(vertex_offset_pattern, s),
            ];
    return select_columns(header, data, matcher);

def extract_material_parameter(header, data):
    young_pattern = "young";
    poisson_pattern = "poisson";
    shear_pattern = "shear";
    #elasticity_mode_pattern = "elasticity_mode";
    matcher = [
            lambda s: re.search(young_pattern, s),
            lambda s: re.search(poisson_pattern, s),
            lambda s: re.search(shear_pattern, s),
            #lambda s: re.search(elasticity_mode_pattern, s),
            ];
    return select_columns(header, data, matcher);

def linear_coeffs(x):
    return x;

def quadratic_coeffs(x):
    x_dof = x.shape[1];
    coeffs = [];

    for i in range(x_dof):
        for j in range(x_dof):
            i_col = x[:, i];
            j_col = x[:, j];
            coeffs.append(np.multiply(i_col, j_col).reshape((-1, 1)));
    return np.hstack(coeffs);

def poly_fit(x, y, degree):
    n = x.shape[0];
    assert(n == y.shape[0]);

    x_dof = x.shape[1];
    y_dof = y.shape[1];

    num_sampled = int(n * 0.025);
    print("{} samples used".format(num_sampled));
    sample_indices = np.random.random_integers(0, n, num_sampled);
    sampled_x = x[sample_indices];
    sampled_y = y[sample_indices];

    coeff_mat = [np.ones((num_sampled, 1))];
    coeff_mat.append(linear_coeffs(sampled_x));
    coeff_mat.append(quadratic_coeffs(sampled_x));

    coeff_mat = np.hstack(coeff_mat);
    sol, residual, rank, singular_vals = lstsq(coeff_mat, sampled_y);

    full_coeff_mat = [np.ones((n, 1))];
    full_coeff_mat.append(linear_coeffs(x));
    full_coeff_mat.append(quadratic_coeffs(x));
    full_coeff_mat = np.hstack(full_coeff_mat);

    err = np.absolute(full_coeff_mat.dot(sol) - y);
    return sol, err;

def generate_function_signature(num_dof, c, values, x_idx):
    linear_terms = [];
    for i in range(num_dof):
        if i == x_idx:
            linear_terms.append("x * {}".format(c[i+1]));
        else:
            linear_terms.append("{}".format(c[i+1] * values[i]));
    linear_terms = "+".join(linear_terms);

    quadratic_terms = [];
    for i in range(num_dof):
        if i == x_idx:
            xi_val = "x";
        else:
            xi_val = values[i];
        for j in range(num_dof):
            if j == x_idx:
                xj_val = "x";
            else:
                xj_val = values[j];

            quadratic_terms.append("{} * {} * {}".format(xi_val, xj_val,
                    c[1+num_dof+i*num_dof+j]));
    quadratic_terms = "+".join(quadratic_terms);

    signature = "\nf <- function(x) {{ {} + {} + {}; }}\n".format(
            c[0], linear_terms, quadratic_terms);
    return signature;

def evaluate(values, c):
    num_dof = len(values);
    r = c[0];

    for i in range(num_dof):
        r += c[i+1] * values[i];

    for i in range(num_dof):
        for j in range(num_dof):
            r += c[1+num_dof+i*num_dof+j] * values[i] * values[j];
    return r;

def dense_sample_slice(num_samples, values, min_values, max_values, c, idx_0, idx_1, idx_2):
    samples_0 = np.linspace(min_values[idx_0], max_values[idx_0], num_samples);
    samples_1 = np.linspace(min_values[idx_1], max_values[idx_1], num_samples);
    samples_2 = np.linspace(min_values[idx_2], max_values[idx_2], num_samples);

    values_0, values_1, values_2 = np.meshgrid(samples_0, samples_1, samples_2);
    values_0 = values_0.ravel();
    values_1 = values_1.ravel();
    values_2 = values_2.ravel();

    data = [];
    for v0,v1,v2 in zip(values_0, values_1, values_2):
        values[idx_0] = v0;
        values[idx_1] = v1;
        values[idx_2] = v2;
        r = evaluate(values, c);
        data.append([v0, v1, v2, r]);

    return data;

def save_slice(filename, header, data):
    with open(filename, 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(header);
        for row in data:
            writer.writerow(["{:10.6f}".format(entry) for entry in row]);

def plot_err_histogram(err, material_data, material_header):
    assert(err.shape == material_data.shape);
    rows = material_data.shape[0];
    err = err.ravel(order="F");
    val = material_data.ravel(order="F");
    rel_err = np.divide(err, val);
    header = np.repeat(material_header, rows);

    err_header = ["{}_err".format(name) for name in header];
    assert(len(err_header) == len(rel_err));
    with open("tmp_err.csv", 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(("relative_err", "type"));
        for err, name in zip(rel_err, err_header):
            writer.writerow(["{:10.6f}".format(err), name]);

    r_script = """
library(ggplot2);
data <- read.csv("tmp_err.csv");
p <- ggplot(data);
p <- p + geom_histogram(aes(x=relative_err), binwidth=0.01);
p <- p + facet_wrap(~ type, ncol=4);
ggsave("tmp_err.pdf", width=10, height=6);
"""

    with open("tmp.r", 'w') as fout:
        fout.write(r_script);

    command = "Rscript tmp.r";
    check_call(command.split());

def plot_slice(csv_file, outer_index, inner_index, color_index, output_index,
        x_header, y_header, values, output_name):
    r_script = """
library(ggplot2);
data <- read.csv("{csv_file}");
slice_data <- read.csv("{slice_file}");
""".format(csv_file = csv_file, slice_file = "tmp_slice.csv");

    for i,header in enumerate(x_header):
        if i == outer_index or i == inner_index or i == color_index:
            continue;
        r_script += "\ndata <- subset(data, {}=={});\n".format(header, values[i]);

    r_script += "\ndata <- subset(data, select=c({}, {}, {}, {}))\n".format(
            x_header[outer_index],
            x_header[inner_index],
            x_header[color_index],
            y_header[output_index]);

    r_script += """
p <- ggplot(data);
p <- p + geom_point(aes(x={x_name}, y={y_name}, color={color_name}));
p <- p + geom_line(aes(x={x_name}, y={y_name}, color={color_name}, group={color_name}));
p <- p + geom_line(data=slice_data, aes(x={x_name}, y={y_name}, group={color_name}));
p <- p + facet_wrap(~ {facet_name}, ncol=5);
p <- p + scale_color_gradientn(colours=rainbow(7));
""".format(x_name = x_header[inner_index],
        y_name = y_header[output_index],
        color_name = x_header[color_index],
        facet_name = x_header[outer_index]);

    r_script += "\nggsave(\"{}\", width=15, height=6);\n"\
            .format(output_name);

    with open("tmp.r", 'w') as fout:
        fout.write(r_script);

    command = "Rscript tmp.r";
    check_call(command.split());

def parse_args():
    parser = argparse.ArgumentParser(description=
            "Find the best fit polynomial that maps pattern parameter to material properties");
    parser.add_argument("--degree", "-d", type=int, default=2,
            help="max degree of the polynomial");
    parser.add_argument("index_dir",
            help="index directory that contains fit.csv");
    return parser.parse_args();

def main():
    args = parse_args();

    csv_file = args.index_dir + "/fit.csv";
    header, data = load_csv(csv_file);

    pattern_header, pattern_data = extract_pattern_parameter(header, data);
    material_header, material_data = extract_material_parameter(header, data);

    sol, err = poly_fit(pattern_data, material_data, args.degree);

    plot_err_histogram(err, material_data, material_header);

    in_idx_0 = 0;
    in_idx_1 = 1;
    in_idx_2 = 4;
    out_idx = 5;
    slice_idx = 20;

    fitted_slice = dense_sample_slice(5, pattern_data[slice_idx],
            np.amin(pattern_data, axis=0),
            np.amax(pattern_data, axis=0),
            sol.T[out_idx], in_idx_0, in_idx_1, in_idx_2);

    save_slice("tmp_slice.csv", [
        pattern_header[in_idx_0],
        pattern_header[in_idx_1],
        pattern_header[in_idx_2],
        material_header[out_idx]],
        fitted_slice);

    plot_slice(csv_file, in_idx_0, in_idx_1, in_idx_2, out_idx,
            pattern_header, material_header, pattern_data[slice_idx], "tmp.pdf");

if __name__ == "__main__":
    main();
