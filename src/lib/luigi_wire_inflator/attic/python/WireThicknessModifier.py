import json
import numpy as np

from WireModifier import WireModifier

import PyWireInflator2DSetting
import PyWireInflator2D

class WireThicknessModifier(WireModifier):
    def __init__(self, config):
        """ Syntax:
        {
            "type": "vertex_orbit" or "edge_orbit",
            "effective_orbits": [i0, i1, ...],
            "thickness": [#, #, ..., "{x} + {y}", ...],
            "default": float 
        }
        """
        self.config = config;
        self.thickness_type = config["type"];
        assert(self.thickness_type == "vertex_orbit");
        self.default_thickness = config["default"];
        self.thickness = config["thickness"];

    def modify(self, param, inflator, **kwargs):
        assert(len(param) == inflator.get_num_parameters());
        self.__load_orbits(inflator);
        self.__compute_thickness_orbit_index_map(inflator);
        num_params = len(param);

        effective_thicknesses = self.__compute_thickness(**kwargs);
        for orbit, thickness in zip(self.effective_orbits, effective_thicknesses):
            # Look for the param index.
            orbit = tuple(sorted(orbit));
            param_index = self.orbit_index_map.get(orbit, -1);
            if param_index < 0:
                raise RuntimeError("Vertex orbit {} is not valid".format(
                    orbit));
            param[param_index] = thickness;

        return param;

    def __load_orbits(self, inflator):
        assert(self.thickness_type == "vertex_orbit");
        self.orbits = np.array(inflator.vertex_orbits);
        self.effective_orbits = self.orbits[self.config["effective_orbits"]];

    def __compute_thickness_orbit_index_map(self, inflator):
        target_param_type = PyWireInflator2D.WireInflatorFacade.THICKNESS;
        num_parameters = inflator.get_num_parameters();
        self.orbit_index_map = {};
        for i in range(num_parameters):
            if inflator.get_parameter_type(i) != target_param_type:
                continue;
            vertex_orbit = sorted(inflator.get_affected_vertex_orbit(i).ravel().astype(int));
            self.orbit_index_map[tuple(vertex_orbit)] = i;

    def __compute_thickness(self, **kwargs):
        effective_thickness = [];
        for orbit, thickness in zip(self.effective_orbits, self.thickness):
            if isinstance(thickness, float):
                effective_thickness.append(thickness);
            elif isinstance(thickness, (str, unicode)):
                formula = thickness.format(**kwargs);
                effective_thickness.append(eval(formula));
        return effective_thickness;
