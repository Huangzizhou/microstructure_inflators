import numpy as np

import PyMesh
import PyWires

from core.WireNetwork import WireNetwork
from ParameterHandler import ParameterHandler

class PyWiresTiler(object):
    def set_single_cell(self, vertices, edges):
        wire_network = WireNetwork();
        wire_network.load(vertices, edges);
        self.set_single_cell_from_wire_network(wire_network);

    def set_single_cell_from_wire_network(self, network):
        network.offset(-network.bbox_center);
        self.pattern = network;
        self.pattern_bbox_min = np.amin(self.pattern.vertices, axis=0);
        self.pattern_bbox_max = np.amax(self.pattern.vertices, axis=0);
        self.pattern_bbox_size = self.pattern_bbox_max - self.pattern_bbox_min;

        self.raw_pattern = self.pattern.raw_wires;

        num_vertices = self.pattern.num_vertices;
        for name, value in self.pattern.attributes.iteritems():
            self.raw_pattern.add_attribute(name, len(value) == num_vertices);
            self.raw_pattern.set_attribute(name, value);

    def initialize_parameter_manager(self, parameters):
        param_handler = ParameterHandler(self.pattern);
        self.parameter_manager = param_handler\
                .convert_to_PyWires_parameter_manager(parameters);

    def tile(self, bbox_min, bbox_max, reps, parameters=[]):
        tiler = PyWires.WireTiler(self.raw_pattern);
        if len(parameters) > 0:
            self.initialize_parameter_manager(parameters);
            tiler.with_parameters(self.parameter_manager);
        self.raw_wire_network = tiler.tile_with_guide_bbox(bbox_min, bbox_max, reps);
        self.__apply_vertex_offset();

    def tile_hex_mesh(self, mesh, parameters=[]):
        tiler = PyWires.WireTiler(self.raw_pattern);
        if len(parameters) > 0:
            self.initialize_parameter_manager(parameters);
            tiler.with_parameters(self.parameter_manager);
        self.raw_wire_network = tiler.tile_with_guide_mesh(mesh);
        self.__apply_vertex_offset();

    def __apply_vertex_offset(self):
        vertices = self.raw_wire_network.get_vertices();
        offset = self.raw_wire_network.get_attribute("vertex_offset");
        self.raw_wire_network.set_vertices(vertices + offset);

    @property
    def wire_network(self):
        network = WireNetwork();
        network.load(
                self.raw_wire_network.get_vertices(),
                self.raw_wire_network.get_edges());
        attr_names = self.raw_wire_network.get_attribute_names();
        for name in attr_names:
            value = np.copy(self.raw_wire_network.get_attribute(name));
            if value.shape[1] == 1:
                value = value.ravel();
            if name == "thickness":
                if self.raw_wire_network.is_vertex_attribute(name):
                    name = "vertex_thickness";
                else:
                    name = "edge_thickness";
            network.attributes.add(name, value);
        return network;

