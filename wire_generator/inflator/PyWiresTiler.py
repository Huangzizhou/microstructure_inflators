import numpy as np

import PyMesh
import PyWires

from core.WireNetwork import WireNetwork

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

    def tile(self, bbox_min, bbox_max, reps, parameters):
        tiler = PyWires.WireTiler(self.raw_pattern);
        tiler.with_parameters(parameters.raw_parameters);
        self.raw_wire_network = tiler.tile_with_guide_bbox(bbox_min, bbox_max, reps);
        self.__apply_vertex_offset();

    def tile_with_hex_mesh(self, mesh, parameters):
        tiler = PyWires.WireTiler(self.raw_pattern);
        tiler.with_parameters(parameters.raw_parameters);
        self.raw_wire_network = tiler.tile_with_guide_mesh(mesh);
        self.__apply_vertex_offset();

    def __apply_vertex_offset(self):
        vertices = self.raw_wire_network.get_vertices();
        offset = self.raw_wire_network.get_attribute("vertex_offset");
        self.raw_wire_network.set_vertices(vertices + offset);

    @property
    def wire_network(self):
        tiled_wire_network = WireNetwork();
        tiled_wire_network.load_from_raw(self.raw_wire_network);
        return tiled_wire_network;
