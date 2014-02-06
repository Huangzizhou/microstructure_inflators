#!/usr/bin/env python

import argparse
import json
import os.path

from ggplot import ggplot
from Histogram import Histogram
from ScatterPlot import ScatterPlot

def create_plot(plots, plot_config):
    plot_type = plot_config.get("type", "unknown");
    if plot_type == "histogram":
        plot = Histogram(plots, plot_config);
    elif plot_type == "scatter_plot":
        plot = ScatterPlot(plots, plot_config);
    else:
        raise NotImplementedError("{} plot type is not supported."\
                .format(plot_type));
    return plot;

def plot(config_file):
    with open(config_file, 'r') as fin:
        config = json.load(fin);

    config_dir = os.path.dirname(config_file);

    csv_file = config.get("csv_file");
    if not os.path.isabs(csv_file):
        # Relative file path is relative to the config file location.
        csv_file = os.path.join(config_dir, csv_file);

    plots = ggplot();
    plots.set_data(csv_file);

    for plot_config in config.get("plots"):
        plot_config["config_dir"] = config_dir;
        plot = create_plot(plots, plot_config);
        plot.generate_r_script();

    r_file = os.path.splitext(csv_file)[0] + ".r";
    plots.plot(r_file);

def parse_args():
    parser = argparse.ArgumentParser(description="Plot csv files");
    parser.add_argument("config_file", help="Plot configuration file");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    plot(args.config_file);

if __name__ == "__main__":
    main();

