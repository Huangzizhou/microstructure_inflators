#!/usr/bin/env python
import argparse
import numpy as np
from mako.template import Template
from mako.runtime import Context
import os.path
from subprocess import check_call
from timethis import timethis
from ggplot import ggplot

TMP_DIR = "./";

def update_prefix(prefix, rank, err_bound):
    if len(prefix) > 0:
        prefix += "_"
    prefix = prefix + "p{}_e{}".format(rank, err_bound);
    return prefix;

def plot_error(plot, prefix, rank, err_bound, output_format):
    plot.scatter_plot("Youngs_modulus", "Poisson_ratio", "Error", False);
    plot.title("Rank-{} Laminates (Error < {})\n".format(rank, err_bound) +\
            "Young's modulus vs Poisson ratio");
    #plot.set_range(base_A_young, base_B_young, 0.288,base_A_poisson);
    #plot.set_range(base_A_young, 5, 0.288,base_A_poisson);
    plot.draw_point(base_A_young, base_A_poisson);
    plot.draw_text(base_A_young, base_A_poisson, " A");
    plot.draw_point(base_B_young, base_B_poisson);
    plot.draw_text(base_B_young, base_B_poisson, " B");
    plot.save_plot("{}_young_poisson_error_{}.{}".format(prefix, base_B_young, output_format));

def plot_ratio(plot, prefix, rank, err_bound, output_format):
    ratio_sum = "+".join(["Ratio_{}".format(i) for i in range(rank)]);
    plot.scatter_plot("Youngs_modulus", "Poisson_ratio", ratio_sum, False);
    plot.title("Rank-{} Laminates (Error < {})\n".format(rank, err_bound) +\
            "Young's modulus vs Poisson ratio");
    plot.set_range(base_A_young, base_B_young, 0.288,base_A_poisson);
    plot.draw_point(base_A_young, base_A_poisson);
    plot.draw_text(base_A_young, base_A_poisson, " A");
    plot.draw_point(base_B_young, base_B_poisson);
    plot.draw_text(base_B_young, base_B_poisson, " B");
    plot.save_plot("{}_young_poisson_ratio.{}".format(prefix, output_format));

def plot_angle_diff(plot, prefix, rank, err_bound, output_format):
    plot.add_column("Angle_Diff",\
            "pmin(abs(Angle_0 - Angle_1), pi - abs(Angle_0 - Angle_1))");
    plot.scatter_plot("Youngs_modulus", "Poisson_ratio", "Angle_Diff", False);
    plot.title("Rank-{} Laminates (Error < {})\n".format(rank, err_bound) +\
            "Young's modulus vs Poisson ratio");
    plot.draw_point(base_A_young, base_A_poisson);
    plot.draw_text(base_A_young, base_A_poisson, " A");
    plot.draw_point(base_B_young, base_B_poisson);
    plot.draw_text(base_B_young, base_B_poisson, " B");
    plot.save_plot("{}_young_poisson_angle.{}".format(prefix, output_format));

@timethis
def plot(csv_file, prefix, rank, err_bound, output_format):
    prefix = update_prefix(prefix, rank, err_bound);
    plot = ggplot();
    plot.set_data(csv_file);
    plot.filter_data("Error", 0.0, err_bound);
    #plot.discretize_column("Angle_0");
    #plot.discretize_column("Ratio_0");

    base_A_young = float(plot.get_col("Young_A")[0]);
    base_B_young = float(plot.get_col("Young_B")[0]);
    base_A_poisson = float(plot.get_col("Poisson_A")[0]);
    base_B_poisson = float(plot.get_col("Poisson_B")[0]);

    #plot.scatter_plot("Lambda", "Mu", "Error");
    #plot.title("Rank-{} Laminates (Error < {})\n".format(rank, err_bound) +\
    #        "Lambda vs Mu");
    #plot.save_plot("{}_lambda_mu.{}".format(prefix, output_format));

    #plot.scatter_plot("Youngs_modulus", "Poisson_ratio", "Error", False);
    #plot.title("Rank-{} Laminates (Error < {})\n".format(rank, err_bound) +\
    #        "Young's modulus vs Poisson ratio");
    ##plot.set_range(base_A_young, base_B_young, 0.288,base_A_poisson);
    ##plot.set_range(base_A_young, 5, 0.288,base_A_poisson);
    #plot.draw_point(base_A_young, base_A_poisson);
    #plot.draw_text(base_A_young, base_A_poisson, " A");
    #plot.draw_point(base_B_young, base_B_poisson);
    #plot.draw_text(base_B_young, base_B_poisson, " B");
    ##plot.save_plot("{}_young_poisson_error.{}".format(prefix, output_format,
    ##    width=5, height=3));
    #plot.save_plot("{}_young_poisson_error_{}.{}".format(prefix, base_B_young, output_format));

    #ratio_sum = "+".join(["Ratio_{}".format(i) for i in range(rank)]);
    #plot.scatter_plot("Youngs_modulus", "Poisson_ratio", ratio_sum, False);
    #plot.title("Rank-{} Laminates (Error < {})\n".format(rank, err_bound) +\
    #        "Young's modulus vs Poisson ratio");
    #plot.set_range(base_A_young, base_B_young, 0.288,base_A_poisson);
    #plot.draw_point(base_A_young, base_A_poisson);
    #plot.draw_text(base_A_young, base_A_poisson, " A");
    #plot.draw_point(base_B_young, base_B_poisson);
    #plot.draw_text(base_B_young, base_B_poisson, " B");
    #plot.save_plot("{}_young_poisson_ratio.{}".format(prefix, output_format));

    plot.add_column("Angle_Diff",\
            "pmin(abs(Angle_0 - Angle_1), pi - abs(Angle_0 - Angle_1))");
    plot.scatter_plot("Youngs_modulus", "Poisson_ratio", "Angle_Diff", False);
    plot.title("Rank-{} Laminates (Error < {})\n".format(rank, err_bound) +\
            "Young's modulus vs Poisson ratio");
    plot.draw_point(base_A_young, base_A_poisson);
    plot.draw_text(base_A_young, base_A_poisson, " A");
    plot.draw_point(base_B_young, base_B_poisson);
    plot.draw_text(base_B_young, base_B_poisson, " B");
    plot.save_plot("{}_young_poisson_angle.{}".format(prefix, output_format));

    plot.plot();

def parse_args():
    parser = argparse.ArgumentParser(description="Plot csv files");
    parser.add_argument("csv_file", help="target csv file");
    parser.add_argument("--rank", help="number of laminations",
            type=int, required=True);
    parser.add_argument("--prefix", help="prefix of output file",
            default="");
    parser.add_argument("-E", "--error-bound",
            help="only plot data within this error bound", type=float,
            default=1.0);
    parser.add_argument("--format", default="png");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    plot(args.csv_file, args.prefix, args.rank, args.error_bound, args.format);
    timethis.summarize();

if __name__ == "__main__":
    main();
