import json
import os.path
import numpy as np

import microstructures_setting
from wire_generator.utils.find_file import find_file
from wire_generator.core.WireNetwork import WireNetwork
from wire_generator.parameter.PyParameters import PyParameters

class DofParameter:
    __known_wires = {};
    __known_params = {};

    def __init__(self, config_file):
        self.__load_config(config_file);
        self.__load_wire();
        self.__load_parameters();

    def __load_config(self, config_file):
        self.config_dir = os.path.dirname(config_file);
        with open(config_file, 'r') as fin:
            self.config = json.load(fin);

    def __load_wire(self):
        wire_file = find_file(str(self.config["wire_network"]), self.config_dir);
        if wire_file in DofParameter.__known_wires:
            self.wire_network = DofParameter.__known_wires[wire_file];
        else:
            self.wire_network = WireNetwork();
            self.wire_network.load_from_file(wire_file);
            DofParameter.__known_wires[wire_file] = self.wire_network;

        basename = os.path.basename(wire_file);
        name, ext = os.path.splitext(basename);
        self.pattern_name = name;

    def __load_parameters(self):
        assert("dof_file" in self.config);
        dof_file = find_file(str(self.config["dof_file"]), self.config_dir);
        if self.pattern_name in DofParameter.__known_params:
            self.parameters = DofParameter.__known_params[self.pattern_name];
            self.parameters.raw_parameters.load_dofs(dof_file);
        else:
            print("param not found");
            self.parameters = PyParameters(
                    self.wire_network, self.config["thickness"]);
            self.parameters.load_dof_file(dof_file);
            DofParameter.__known_params[self.pattern_name] = self.parameters;
        self.dofs = np.copy(self.parameters.dofs);


