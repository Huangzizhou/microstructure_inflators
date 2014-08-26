import json
import numpy as np
from numpy.linalg import norm

from WireModifier import WireModifier

class WireVertexOffsetModifier(WireModifier):
    def __init__(self, config):
        """ Syntax:
        {
            "type": "vertex_orbit",
            "orbit_file": "filename",
            "effective_orbits": [i0, i1, ...],
            "offset_percentages": [[#, #, #], [0.0, #, 0.0], ...,\
                    ["{x}", "{y}", 0.0]]
        }
        """
        self.offset_type = config["type"];
        self.effective_orbits = config["effective_orbits"];
        self.offset_percentages = config["offset_percentages"];

    def modify(self, param, **kwargs):
        offsets = np.zeros(5);
        self.__compute_offset(offsets, **kwargs);
        param[5:] = np.array(offsets)[[0,2,4]];
        return param;

    def __compute_offset(self, offsets, **kwargs):
        for orbit, offset_percent in zip(
                self.effective_orbits, self.offset_percentages):
            offset = self.__compute_offset_per_orbit(orbit, offset_percent, **kwargs);
            offsets[orbit] = offset;

    def __compute_offset_per_orbit(
            self, orbit, offset_percent, **kwargs):
        offset_percent = [eval(p.format(**kwargs)) if isinstance(p, str) else p
                for p in offset_percent];
        return norm(offset_percent);
