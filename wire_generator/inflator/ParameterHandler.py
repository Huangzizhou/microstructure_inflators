import numpy as np
from parameter.VertexThicknessParameter import VertexThicknessParameter
from parameter.EdgeThicknessParameter import EdgeThicknessParameter
from parameter.VertexOffsetParameter import VertexOffsetParameter

class ParameterHandler(object):
    def __init__(self, wire_network):
        self.wire_network = wire_network;
        self.__remove_old_attributes();
        self.__initialize_orbits();

    def __remove_old_attributes(self):
        if "vertex_thickness" in self.wire_network.attributes:
            del self.wire_network.attributes["vertex_thickness"];
        if "edge_thickness" in self.wire_network.attributes:
            del self.wire_network.attributes["edge_thickness"];
        if "vertex_offset" in self.wire_network.attributes:
            del self.wire_network.attributes["vertex_offset"];

    def __initialize_orbits(self):
        if "symmetry_vertex_orbit" not in self.wire_network.attributes or\
                "symmetry_edge_orbit" not in self.wire_network.attributes:
            self.wire_network.compute_symmetry_orbits();

        self.vertex_orbits = self.wire_network.attributes["symmetry_vertex_orbit"];
        self.edge_orbits = self.wire_network.attributes["symmetry_edge_orbit"];

    def convert_to_attributes(self, parameters, **kwargs):
        attributes = self.wire_network.attributes;
        for param in parameters:
            if isinstance(param, VertexThicknessParameter):
                if "vertex_thickness" not in attributes:
                    attributes["vertex_thickness"] = np.zeros(
                            self.wire_network.num_vertices);
                thickness = attributes["vertex_thickness"];
                affected_vertices = self.vertex_orbits == param.orbit_id;
                thickness[affected_vertices] = param.evaluate(**kwargs);
            elif isinstance(param, EdgeThicknessParameter):
                if "edge_thickness" not in attributes:
                    attributes["edge_thickness"] = np.zeros(
                            self.wire_network.num_edges);
                thickness = attributes["edge_thickness"];
                affected_edges = self.edge_orbits == param.orbit_id;
                thickness[affected_edges] = param.evaluate(**kwargs);
            elif isinstance(param, VertexOffsetParameter):
                if "vertex_offset" not in attributes:
                    attributes["vertex_offset"] = np.zeros((
                        self.wire_network.num_vertices,
                        self.wire_network.dim));
                offset = attributes["vertex_offset"];
                affected_vertices = self.vertex_orbits == param.orbit_id;
                percentage = param.evaluate(**kwargs);
                vertices = self.wire_network.vertices[affected_vertices];
                centroid = np.mean(vertices, axis=0);
                offset[affected_vertices] = (vertices - centroid) * percentage;
            else:
                raise NotImplementedError(
                        "Unknown parameter type: {}".format(param));
