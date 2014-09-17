from WireModifier import WireModifier

import json
import numpy as np

class WireThicknessModifier(WireModifier):
    def __init__(self, config):
        """ Syntax:
        "thickness": {
            "type": "vertex_orbit" or "edge_orbit",
            "effective_orbits": [i0, i1, ...],
            "thickness": [#, #, ..., "{x} + {y}", ...],
            "default": float 
        }
        """
        self.config = config;
        self.thickness_type = config["type"];
        self.thicknesses = config["thickness"];
        self.default_thickness = config["default"];

    def modify(self, wire_network, **kwargs):
        self.__load_orbits(wire_network);
        thicknesses = self.__generate_default_thicknesses(wire_network);
        self.__compute_thickness(wire_network, thicknesses, **kwargs);
        self.__assign_thickness_attribute(wire_network, thicknesses);

    def __load_orbits(self, wire_network):
        if "symmetry_vertex_orbit" not in wire_network.attributes or\
                "symmetry_edge_orbit" not in wire_network.attributes:
            wire_network.compute_symmetry_orbits();
        if self.thickness_type == "vertex_orbit":
            self.orbits = wire_network.attributes["symmetry_vertex_orbit"];
        elif self.thickness_type == "edge_orbit":
            self.orbits = wire_network.attributes["symmetry_edge_orbit"];
        else:
            raise NotImplementedError("Thickness type ({}) is not supported"\
                    .format(self.thickness_type));

        self.orbits = np.array(self.orbits);
        self.effective_orbits = self.orbits[self.config["effective_orbits"]];

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

