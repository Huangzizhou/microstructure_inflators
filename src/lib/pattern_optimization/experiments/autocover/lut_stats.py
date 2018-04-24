#!/usr/bin/env python
# Compute statistics of lookup table values.
# (Could be used on individual Lipschitz components).
import sys
from LookupTable import LUT

t = LUT(sys.argv[1])
if (len(sys.argv) == 3):
    t.filterPattern(int(sys.argv[2]))

t.reportStats()
