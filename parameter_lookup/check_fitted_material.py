#!/usr/bin/env python

import argparse
import csv
import json
import numpy as np
import os.path
from subprocess import check_call

import PyMeshSetting
import PyMesh

def load_mesh(mesh_file):
    factory = PyMesh.MeshFactory();
    factory.load_file(mesh_file);
    factory.drop_zero_dim();
    return factory.create();

def parse_args():
    parser = argparse.ArgumentParser(
            description="Check the distribution of Young's modulus");
    parser.add_argument("input_mesh", help="input mesh");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();

    basename, ext = os.path.splitext(args.input_mesh);
    path, name = os.path.split(basename);

    mesh = load_mesh(args.input_mesh);

    target_young = mesh.get_attribute("young").ravel();
    fitted_young = mesh.get_attribute("fitted_young_x").ravel();

    TMP_DIR = "/tmp/"

    csv_file = os.path.join(TMP_DIR, "{}.csv".format(name));
    with open(csv_file, 'w') as fout:
        writer = csv.writer(fout);
        writer.writerow(["type", "value"]);
        for target, fitted in zip(target_young, fitted_young):
            writer.writerow(["fitted", fitted]);
            writer.writerow(["target", target]);

    png_name = basename + "_young.png";
    r_script = """
    library(ggplot2);
    data <- read.csv("{}");
    data$value = as.numeric(as.character(data$value));
    p <- ggplot(data);
    p <- p + geom_histogram(aes(x=value, y=..density.., fill=type), binwidth=10) + geom_density(aes(x=value));
    p <- p + facet_wrap(~ type);
    p <- p + xlab("X Young's modulus");
    ggsave(file="{}", width=10, height=6);
    """.format(csv_file, png_name);

    tmp_r_file = os.path.join(TMP_DIR, "{}.r".format(name));
    with open(tmp_r_file, 'w') as fout:
        fout.write(r_script);

    command = "Rscript {}".format(tmp_r_file);
    check_call(command.split());

if __name__ == "__main__":
    main();
