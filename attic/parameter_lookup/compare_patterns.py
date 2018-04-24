#!/usr/bin/env python

import argparse
import re
import os
import os.path
from subprocess import check_call

load_csv_script = """
data_{pattern_name} <- read.csv("{csv_file}");
data_{pattern_name} <- subset(data_{pattern_name},
    select=c(young_x, young_y, young_z, shear_yz,
             shear_zx, shear_xy, poisson_yz,
             poisson_zx, poisson_xy));
data_{pattern_name}$pattern_name <- "{pattern_name}";
"""

r_script = """
library(ggplot2);
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
"#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
{load_script}
data <- rbind({data_frames});

p <- ggplot(data);
p <- p + geom_point(aes(young_x, poisson_xy, color=shear_xy));
p <- p + scale_colour_gradientn( colours = jet.colors(7));
p <- p + facet_wrap(~pattern_name, ncol={num_patterns});
p <- p + ggtitle("X Young's modulus vs XY Poisson ratio");
ggsave("{output_file}", width=4*{num_patterns}, height=5);
"""

def parse_args():
    parser = argparse.ArgumentParser(
            description="Compare homogenized results across different patterns");
    parser.add_argument("--output", "-o", help="output file");
    parser.add_argument("index_dir", nargs="+",
            help="index directory of patterns");
    return parser.parse_args();

def main():
    args = parse_args();

    index_pattern = re.compile("\dD/(.*)/.*mm_cell/");
    pattern_names = [];
    load_scripts = [];
    for index_dir in args.index_dir:
        csv_file = os.path.join(index_dir, "fit.csv");
        assert(os.path.exists(csv_file));
        m = index_pattern.search(os.path.abspath(csv_file));
        assert(m is not None);
        pattern_name = m.group(1);
        print("processing {}".format(pattern_name));

        pattern_names.append(pattern_name);
        load_scripts.append(load_csv_script.format(
            pattern_name = pattern_name,
            csv_file = csv_file));

    script = r_script.format(
            load_script="\n".join(load_scripts),
            data_frames=",".join(["data_{}".format(name) for name in pattern_names]), 
            num_patterns=len(set(pattern_names)),
            output_file=args.output);

    r_file = "tmp.r";
    with open(r_file, 'w') as fout:
        fout.write(script);

    command = "Rscript {}".format(r_file);
    check_call(command.split());

    os.remove(r_file);

if __name__ == "__main__":
    main();
