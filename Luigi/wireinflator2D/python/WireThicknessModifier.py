import json
import numpy as np

from WireModifier import WireModifier

class WireThicknessModifier(WireModifier):
    def __init__(self, config):
        """ Syntax:
        {
            "type": "vertex_orbit" or "edge_orbit",
            "orbit_file": "filename",
            "effective_orbits": [i0, i1, ...],
            "thickness": [#, #, ..., "{x} + {y}", ...],
            "default": float 
        }
        """
        self.thickness_type = config["type"];
        assert(self.thickness_type == "vertex_orbit");
        self.effective_orbits = config["effective_orbits"];
        self.default_thickness = config["default"];
        self.thickness = config["thickness"];

    def modify(self, param, **kwargs):
        thickness = np.ones(5) * self.default_thickness;
        self.__compute_thickness(thickness, **kwargs);
        param[:5] = thickness;

    def __compute_thickness(self, thicknesses, **kwargs):
        for orbit, thickness in zip(self.effective_orbits, self.thickness):
            if isinstance(thickness, float):
                thicknesses[orbit] = thickness;
            elif isinstance(thickness, (str, unicode)):
                formula = thickness.format(**kwargs);
                thicknesses[orbit] = eval(formula);
