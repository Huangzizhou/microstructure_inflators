import numpy as np
from numpy.linalg import norm

from WireReader import WireReader

class WireNetwork(object):

    def load(self, vertices, edges):
        self.vertices = np.array(vertices);
        self.edges = np.array(edges);
        self.dim = self.vertices.shape[1];
        self.__compute_connectivity();

    def load_from_file(self, wire_file):
        self.__parse_wire_file(wire_file);
        self.__compute_connectivity();

    def scale(self, factors):
        if isinstance(factors, float):
            self.vertices *= factors;
        else:
            self.vertices = np.multiply(self.vertices, factors);

    def translate(self, offset):
        self.vertices += offset;

    def trim(self):
        """ Remove all hanging edges
        e.g. edge with at least one vertex of valance <= 1
        """
        edge_to_keep = np.all(self.vertex_valance[self.edges] > 1, axis=1);
        self.edges = self.edges[edge_to_keep];
        self.__remove_isolated_vertices();
        self.__compute_connectivity();
        if len(self.vertices) == 0:
            raise RuntimeError("Zero vertices left after trimming.");

    def __parse_wire_file(self, wire_file):
        parser = WireReader(wire_file);
        self.dim = parser.dim;
        self.vertices = np.array(parser.vertices);
        self.edges = np.array(parser.edges);

    def __compute_connectivity(self):
        self.vertex_neighbors = [[] for i in range(self.num_vertices)];
        self.vertex_edge_neighbors = [[] for i in range(self.num_vertices)];
        self.vertex_valance = np.zeros(self.num_vertices);
        for i,edge in enumerate(self.edges):
            self.vertex_neighbors[edge[0]].append(edge[1]);
            self.vertex_neighbors[edge[1]].append(edge[0]);
            self.vertex_edge_neighbors[edge[0]].append(i);
            self.vertex_edge_neighbors[edge[1]].append(i);
            self.vertex_valance[edge[0]] += 1;
            self.vertex_valance[edge[1]] += 1;

    def __remove_isolated_vertices(self):
        vertex_map = np.zeros(self.num_vertices, dtype=int) - 1;
        vertices = [];
        for edge in self.edges:
            if vertex_map[edge[0]] == -1:
                vertex_map[edge[0]] = len(vertices);
                vertices.append(self.vertices[edge[0]]);
            if vertex_map[edge[1]] == -1:
                vertex_map[edge[1]] = len(vertices);
                vertices.append(self.vertices[edge[1]]);

        self.vertices = np.array(vertices);
        self.edges = vertex_map[self.edges];

    @property
    def num_vertices(self):
        return len(self.vertices);

    @property
    def num_edges(self):
        return len(self.edges);

    @property
    def bbox(self):
        return np.amin(self.vertices, axis=0), np.amax(self.vertices, axis=0);

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

