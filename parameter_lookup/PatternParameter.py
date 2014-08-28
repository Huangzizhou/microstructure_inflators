import json
import os.path

import numpy as np
from numpy.linalg import norm

class PatternParameter:
    def __init__(self, dim, modifier_file):
        self.dim = dim;
        self.root_dir = os.path.dirname(modifier_file);

        with open(modifier_file, 'r') as fin:
            config = json.load(fin);
            config["root_dir"] = self.root_dir;
            if "thickness" in config:
                self.thickness_parameter = \
                        ThicknessParameter(config);
            if "vertex_offset" in config:
                self.vertex_offset_parameter =\
                        VertexOffsetParameter(dim, config);

    @property
    def names(self):
        names = [];
        if hasattr(self, "thickness_parameter"):
            names += self.thickness_parameter.names;
        if hasattr(self, "vertex_offset_parameter"):
            names += self.vertex_offset_parameter.names;
        return names;

    @property
    def values(self):
        values = [];
        if hasattr(self, "thickness_parameter"):
            values += self.thickness_parameter.values;
        if hasattr(self, "vertex_offset_parameter"):
            values += self.vertex_offset_parameter.values;
        return values;

class ThicknessParameter(PatternParameter):
    def __init__(self, config):
        self.root_dir = config["root_dir"]
        self.thickness_config = config["thickness"];
        self.load_orbits();
        self.load_thickness_parameters();

    def load_orbits(self):
        orbit_file = os.path.join(self.root_dir,
                self.thickness_config["orbit_file"]);
        with open(orbit_file, 'r') as fin:
            self.orbit_config = json.load(fin);
            self.num_vertex_orbits = len(self.orbit_config["vertex_orbits"]);
            self.num_edge_orbits = len(self.orbit_config["edge_orbits"]);

    def load_thickness_parameters(self):
        if self.thickness_config["type"] == "edge_orbit":
            self.load_edge_thickness_parameters();
        else:
            self.load_vertex_thickness_parameters();

    def load_edge_thickness_parameters(self):
        self.values = np.ones(self.num_edge_orbits) *\
                self.thickness_config["default"];
        effective_orbits = self.thickness_config["effective_orbits"];
        thickness = self.thickness_config["thickness"];
        self.values[effective_orbits] = thickness;

        self.values = self.values.tolist();
        self.names = ["edge_orbit_{}_thickness".format(i) for i in
                range(self.num_edge_orbits)];

    def load_vertex_thickness_parameters(self):
        self.values = np.ones(self.num_vertex_orbits) *\
                self.thickness_config["default"];
        effective_orbits = self.thickness_config["effective_orbits"];
        thickness = self.thickness_config["thickness"];
        self.values[effective_orbits] = thickness;

        self.values = self.values.tolist();
        self.names = ["vertex_orbit_{}_thickness".format(i) for i in
                range(self.num_vertex_orbits)];

class VertexOffsetParameter(PatternParameter):
    def __init__(self, dim, config):
        self.dim = dim;
        self.root_dir = config["root_dir"];
        self.offset_config = config["vertex_offset"];
        self.load_orbits();
        self.load_offset_parameters();

    def load_orbits(self):
        orbit_file = os.path.join(self.root_dir,
                self.offset_config["orbit_file"]);
        with open(orbit_file, 'r') as fin:
            self.orbit_config = json.load(fin);
            self.num_vertex_orbits = len(self.orbit_config["vertex_orbits"]);
            self.num_edge_orbits = len(self.orbit_config["edge_orbits"]);

    def load_offset_parameters(self):
        self.values = np.zeros((self.num_vertex_orbits, self.dim));
        effective_orbits = self.offset_config["effective_orbits"];
        offsets = self.offset_config["offset_percentages"];

        self.values[effective_orbits] = offsets;
        self.values = self.values.ravel(order="C").tolist();
        self.names = ["vertex_orbit_{}_offset_{}".format(i, j)
                for j in range(self.dim)
                for i in range(self.num_vertex_orbits) ];

