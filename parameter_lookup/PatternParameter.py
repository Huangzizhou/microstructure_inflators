import json
import os.path

import numpy as np
from numpy.linalg import norm

import microstructures_setting
from wire_generator.utils.find_file import find_file
from wire_generator.core.WireNetwork import WireNetwork
from wire_generator.parameter.ParameterFactory import ParameterFactory
from wire_generator.parameter.EdgeThicknessParameter import EdgeThicknessParameter
from wire_generator.parameter.VertexOffsetParameter import VertexOffsetParameter
from wire_generator.parameter.VertexThicknessParameter import VertexThicknessParameter

class PatternParameter:
    def __init__(self, config_file):
        self.__load_config(config_file);
        self.__load_wire();
        self.__load_parameters();

    def __load_config(self, config_file):
        self.config_dir = os.path.dirname(config_file);
        with open(config_file, 'r') as fin:
            self.config = json.load(fin);

    def __load_wire(self):
        wire_file = find_file(self.config["wire_network"], self.config_dir);
        self.wire_network = WireNetwork();
        self.wire_network.load_from_file(wire_file);

    def __load_parameters(self):
        modifier_file = self.config.get("modifier_file");
        if modifier_file is not None:
            modifier_file = find_file(modifier_file, self.config_dir);
        factory = ParameterFactory(self.wire_network, self.config["thickness"]);
        factory.create_parameters_from_file(modifier_file);
        self.parameters = factory.parameters;

    @property
    def names(self):
        names = [];
        for param in self.parameters:
            names += param.names;
        return names;

    @property
    def values(self):
        values = [];
        for param in self.parameters:
            values += param.evaluate().tolist();
        return values;

    @property
    def modifier_config(self):
        thickness_config = {
                "type": "vertex_orbit",
                "effective_orbits": [],
                "thickness": [],
                "default": 0.5
                };
        offset_config ={
                    "type": "vertex_orbit",
                    "effective_orbits": [],
                    "offset_percentages": []
                };
        for param in self.parameters:
            if isinstance(param, VertexThicknessParameter):
                thickness_config["type"] = "vertex_orbit"
                thickness_config["effective_orbits"].append(param.orbit_id);
                thickness_config["thickness"].append("{{{}}}".format(param.names[0]));
            elif isinstance(param, EdgeThicknessParameter):
                thickness_config["type"] = "edge_orbit"
                thickness_config["effective_orbits"].append(param.orbit_id);
                thickness_config["thickness"].append("{{{}}}".format(param.names[0]));
            elif isinstance(param, VertexOffsetParameter):
                offset_config["effective_orbits"].append(param.orbit_id);
                offset_values = ["{{{}}}".format(val) for val in param.names];
                offset_config["offset_percentages"].append(offset_values);
            else:
                raise NotImplementedError("Unknown parameter of type {}".format(
                    type(param)));

        return {"thickness": thickness_config, "vertex_offset": offset_config};

