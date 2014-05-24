from WireModifier import WireModifier

import json
import numpy as np

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
        self.__load_orbits(config["orbit_file"]);
        self.effective_orbits = self.orbits[config["effective_orbits"]];
        self.thicknesses = config["thickness"];
        self.default_thickness = config["default"];

    def modify(self, wire_network, **kwargs):
        thicknesses = self.__generate_default_thicknesses(wire_network);
        self.__compute_thickness(wire_network, thicknesses, **kwargs);
        self.__assign_thickness_attribute(wire_network, thicknesses);

    def __load_orbits(self, orbit_file):
        with open(orbit_file, 'r') as fin:
            contents = json.load(fin);
            if self.thickness_type == "vertex_orbit":
                orbits = contents["vertex_orbits"];
            elif self.thickness_type == "edge_orbit":
                orbits = contents["edge_orbits"];
            else:
                raise NotImplementedError("Thickness type ({}) is not supported"\
                        .format(self.thickness_type));
            self.orbits = np.array(orbits);

    def __generate_default_thicknesses(self, wire_network):
        if self.thickness_type == "vertex_orbit":
            thicknesses = np.ones(wire_network.num_vertices);
        elif self.thickness_type == "edge_orbit":
            thicknesses = np.ones(wire_network.num_edges);
        else:
            raise NotImplementedError("Thickness type ({}) is not supported"\
                    .format(self.thickness_type));
        return thicknesses * self.default_thickness;

    def __compute_thickness(self, wire_network, thicknesses, **kwargs):
        for orbit, thickness in zip(self.effective_orbits, self.thicknesses):
            if isinstance(thickness, float):
                thicknesses[orbit] = thickness;
            elif isinstance(thickness, (str, unicode)):
                formula = thickness.format(**kwargs);
                thicknesses[orbit] = eval(formula);

    def __assign_thickness_attribute(self, wire_network, thicknesses):
        if "vertex_thickness" in wire_network.attributes:
            del wire_network.attributes["vertex_thickness"];
        if "edge_thickness" in wire_network.attributes:
            del wire_network.attributes["edge_thickness"];

        if self.thickness_type == "vertex_orbit":
            attr_name = "vertex_thickness";
        elif self.thickness_type == "edge_orbit":
            attr_name = "edge_thickness";
        else:
            raise NotImplementedError("Thickness type ({}) is not supported"\
                    .format(self.thickness_type));

        if attr_name not in wire_network.attributes:
            wire_network.attributes.add(attr_name, thicknesses);
        else:
            wire_network.attributes[attr_name] = thicknesses;

