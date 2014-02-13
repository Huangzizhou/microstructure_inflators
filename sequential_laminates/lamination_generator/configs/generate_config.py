#!/usr/bin/env python
import numpy as np
from mako.template import Template

def main():
    config = Template(filename="config.mako");
    for i,area in enumerate(np.arange(0.00001, 0.001, 0.00002)):
        contents = config.render(index=i, max_area=area);
        config_filename = "rank_2_laminate_{}.config".format(i);
        with open(config_filename, 'w') as fout:
            fout.write(contents);

if __name__ == "__main__":
    main();
