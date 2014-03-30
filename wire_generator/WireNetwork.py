import numpy as np

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

    def __parse_wire_file(self, wire_file):
        parser = WireReader(wire_file);
        self.dim = parser.dim;
        self.vertices = np.array(parser.vertices);
        self.edges = np.array(parser.edges);

    def __compute_connectivity(self):
        self.vertex_neighbors = [[] for i in range(self.num_vertices)];
        self.vertex_edge_neighbors = [[] for i in range(self.num_vertices)];
        for i,edge in enumerate(self.edges):
            self.vertex_neighbors[edge[0]].append(edge[1]);
            self.vertex_neighbors[edge[1]].append(edge[0]);
            self.vertex_edge_neighbors[edge[0]].append(i);
            self.vertex_edge_neighbors[edge[1]].append(i);

    @property
    def num_vertices(self):
        return len(self.vertices);

    @property
    def num_edges(self):
        return len(self.edges);

    @property
    def bbox(self):
        return np.amin(self.vertices, axis=0), np.amax(self.vertices, axis=0);
