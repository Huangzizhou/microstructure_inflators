#!/usr/bin/env python

import argparse
import os
import os.path
from subprocess import check_call

r_script = """
library(ggplot2);
data <- read.csv("{csv_file}");

summary(subset(data, select=c(young_x, young_y, young_z)));
summary(subset(data, select=c(shear_yz, shear_zx, shear_xy)));
summary(subset(data, select=c(poisson_yz, poisson_zx, poisson_xy)));

data$pentamode <- data$elasticity_mode_0 / data$elasticity_mode_1;
summary(subset(data, select=pentamode));

p <- ggplot(data);
p <- p + geom_point(aes(young_x, shear_xy));
p <- p + ggtitle("X Young's modulus vs XY shear modulus");
ggsave("{index_dir}/young_x_vs_shear_xy.png", p, width=8, height=6);

p <- ggplot(data);
p <- p + geom_point(aes(young_x, poisson_xy));
p <- p + ggtitle("X Young's modulus vs XY Poisson ratio");
ggsave("{index_dir}/young_x_vs_poisson_xy.png", p, width=8, height=6);

p <- ggplot(data);
p <- p + geom_point(aes(young_x, young_y));
p <- p + ggtitle("X Young's modulus vs Y Young's modulus");
ggsave("{index_dir}/young_x_vs_young_y.png", p, width=8, height=6);
"""

def parse_args():
    parser = argparse.ArgumentParser(
            description="plot material property distribution");
    parser.add_argument("--index-dir", help="index directory");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    csv_file = os.path.join(args.index_dir, "fit.csv");
    contents = r_script.format(csv_file = csv_file, index_dir = args.index_dir);
    r_file = "tmp.r";
    with open(r_file, 'w') as fout:
        fout.write(contents);

    command = "Rscript {}".format(r_file);
    check_call(command.split());

    os.remove(r_file);

if __name__ == "__main__":
    main();
