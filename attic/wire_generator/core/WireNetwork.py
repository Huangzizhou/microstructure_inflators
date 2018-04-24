import numpy as np
from numpy.linalg import norm

import PyMeshSetting
import PyWires

class WireNetwork(object):
    def __init__(self):
        self.raw_wires = PyWires.WireNetwork();
        self.__initialize_wires();

    def load(self, vertices, edges):
        self.raw_wires = PyWires.WireNetwork.create_raw(vertices, edges);
        self.__initialize_wires();

    def load_from_file(self, wire_file):
        self.raw_wires = PyWires.WireNetwork.create(wire_file);
        self.__initialize_wires();

    def load_from_raw(self, raw_wires):
        self.raw_wires = raw_wires;
        self.__initialize_wires();

    def scale(self, factors):
        if isinstance(factors, float):
            factors = np.ones(self.dim) * factors;
        self.raw_wires.scale(factors);

    def offset(self, offset_vector):
        vertices = self.vertices + offset_vector;
        self.vertices = vertices;

    def translate(self, offset):
        self.raw_wires.translate(offset);

    def trim(self):
        """ Remove all hanging edges
        e.g. edge with at least one vertex of valance <= 1
        """
        while np.any(self.vertex_valance <= 1):
            edge_to_keep = np.all(self.vertex_valance[self.edges] > 1,
                    axis=1).tolist();
            self.raw_wires.filter_edges(edge_to_keep);
            vertex_to_keep = [len(self.get_vertex_neighbors(i)) > 0 for i in
                    range(self.num_vertices)];
            self.raw_wires.filter_vertices(vertex_to_keep);

            self.__initialize_wires();
            if len(self.vertices) == 0:
                raise RuntimeError("Zero vertices left after trimming.");

    def compute_symmetry_orbits(self):
        self.add_attribute("vertex_symmetry_orbit");
        self.add_attribute("vertex_cubic_symmetry_orbit");
        self.add_attribute("edge_symmetry_orbit");
        self.add_attribute("edge_cubic_symmetry_orbit");

    def has_attribute(self, name):
        return self.raw_wires.has_attribute(name);

    def add_attribute(self, name, value=None):
        if not self.has_attribute(name):
            self.raw_wires.add_attribute(name);
        if value is not None:
            self.raw_wires.set_attribute(name, value);

    def get_attribute(self, name):
        assert(self.has_attribute(name));
        return self.raw_wires.get_attribute(name);

    def is_vertex_attribute(self, name):
        assert(self.has_attribute(name));
        return self.raw_wires.is_vertex_attribute(name);

    def set_attribute(self, name, value):
        assert(self.has_attribute(name));
        self.raw_wires.set_attribute(name, value);

    def get_vertex_neighbors(self, i):
        if not self.raw_wires.with_connectivity():
            self.raw_wires.compute_connectivity();
        return self.raw_wires.get_vertex_neighbors(i);

    def __initialize_wires(self):
        self.raw_wires.compute_connectivity();
        self.vertex_valance = np.array([
                len(self.raw_wires.get_vertex_neighbors(i)) for i in
                range(self.num_vertices) ], dtype=int);

    @property
    def dim(self):
        return self.raw_wires.get_dim();

    @property
    def num_vertices(self):
        return self.raw_wires.get_num_vertices();

    @property
    def num_edges(self):
        return self.raw_wires.get_num_edges();

    @property
    def vertices(self):
        return self.raw_wires.get_vertices();

    @vertices.setter
    def vertices(self, value):
        self.raw_wires.set_vertices(value);

    @property
    def edges(self):
        return self.raw_wires.get_edges();

    @edges.setter
    def edges(self, value):
        self.raw_wires.set_edges(value);
        self.__initialize_wires();

    @property
    def bbox(self):
        return (self.raw_wires.get_bbox_min().ravel(),
                self.raw_wires.get_bbox_max().ravel());

    @property
    def bbox_center(self):
        return self.raw_wires.center().ravel();

    @property
    def centroid(self):
        return np.average(self.vertices, axis=0);

    @property
    def total_wire_length(self):
        return np.sum(self.wire_lengths);

    @property
    def wire_lengths(self):
        return norm(
            self.vertices[self.edges[:,0]] -
            self.vertices[self.edges[:,1]], axis=1);

    @property
    def attribute_names(self):
        return self.raw_wires.get_attribute_names();

