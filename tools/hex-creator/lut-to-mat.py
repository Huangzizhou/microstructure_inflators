#!/usr/bin/env python
import csv
import sys
import math

import numpy as np
import re

import hexlib

if len(sys.argv) != 3:
    print "usage: ./lut2mat.py <lut table> <matlab table>"
    print "example: ./lut2mat.py sweep.txt sweep.mat"
    sys.exit(-1)


table_path = sys.argv[1]
data_path = sys.argv[2]

out_table = []

with open(table_path, 'r') as table_file:
    csv_reader = csv.reader(table_file, delimiter=' ')


    for row in csv_reader:

        # parse parameters
        m = re.search('#hexa', row[5])
        if m:
            # old type of result file:
            E = float(row[1])
            nu = float(row[2])
            alpha = float(row[3])

            # parse parameters
            m = re.search('-n(.+?)-s', row[5])
            if m:
                n = int(m.group(1))
            else:
                n = 0

            m = re.search('-s(.+?)-t', row[5])
            if m:
                s = float(m.group(1))
            else:
                s = 0.0

            m = re.search('-t(.+?).msh', row[5])
            if m:
                t = float(m.group(1))
            else:
                t = 0.0

            p1 = s/2.0
            p2 = n
            p3 = 1.0
            p4 = (n*t) / s

            volume, volfrac = hexlib.compute_volume_info('positive', p1, p2, p3, p4)


        else:
            # new type of result file
            E = float(row[1])
            nu = float(row[2])
            alpha = float(row[3])

            # parse parameters
            m = re.search('volfrac-(.+?)_', row[5])
            if m:
                volfrac = float(m.group(1))
            else:
                volfrac = 0.0

            m = re.search('p1-(.+?)_', row[5])
            if m:
                p1 = float(m.group(1))
            else:
                p1 = 0.0

            m = re.search('p2-(.+?)_', row[5])
            if m:
                p2 = float(m.group(1))
            else:
                p2 = 0

            m = re.search('p3-(.+?)_', row[5])
            if m:
                p3 = float(m.group(1))
            else:
                p3 = 0.0

            m = re.search('p4-(.+?)_', row[5])
            if m:
                p4 = float(m.group(1))
            else:
                p4 = 0.0

        out_table.append([E, nu, volfrac, p1, p2, p3, p4, alpha])

np.savetxt(data_path, out_table)