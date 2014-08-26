import json
import numpy as np

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
        offsets = np.zeros(3);
        self.__compute_offset(offsets, **kwargs);
        param[5:] = offsets;
        return param;

    def __compute_offset(self, offsets, **kwargs):
        for orbit, offset_percent in zip(
                self.effective_orbits, self.offset_percentages):
            if isinstance(offset_percent, float):
                offsets[orbit] = offset_percent;
            elif isinstance(offset_percent, (str, unicode)):
                offsets[orbit] = eval(offset_percent.format(**kwargs));
