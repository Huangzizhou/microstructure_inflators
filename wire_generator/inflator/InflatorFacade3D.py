import numpy as np
from InflatorFacade import InflatorFacade
from WireTiler import WireTiler
from WireInflator import WireInflator
from PeriodicWireInflator import PeriodicWireInflator
from PyWiresInflator import PyWiresInflator

class InflatorFacade3D(InflatorFacade):
    def __init__(self, wire_network, parameters):
        super(InflatorFacade3D, self).__init__(wire_network, parameters);

    def inflate_with_guide_box(self, bbox_min, bbox_max, repetitions, options):
        tiled_network = self.__tile_with_guide_box(bbox_min, bbox_max,
                repetitions);
        mesh = self.__inflate(tiled_network, options);
        return mesh;

    def inflate_with_guide_mesh(self, mesh, options):
        tiled_network = self.__tile_with_mesh(mesh);
        mesh = self.__inflate(tiled_network, options);
        return mesh;

    def __tile_with_guide_box(self, bbox_min, bbox_max, repetitions):
        bbox_min = np.array(bbox_min);
        bbox_max = np.array(bbox_max);
        unit_bbox_min, unit_bbox_max = self.unit_pattern.bbox;
        unit_bbox_size = unit_bbox_max - unit_bbox_min;
        target_bbox_size = np.divide(bbox_max - bbox_min, repetitions);
        non_zero_dim = unit_bbox_size > 0.0;
        factor = np.ones(self.unit_pattern.dim);
        factor[non_zero_dim] = np.divide(
                target_bbox_size[non_zero_dim],
                unit_bbox_size[non_zero_dim]);
        self.unit_pattern.scale(factor);

        tiler = WireTiler();
        tiler.set_single_cell_from_wire_network(self.unit_pattern);
        tiler.tile(repetitions, self.parameters);

        tiled_network = tiler.wire_network;

        tiled_bbox_min, tiled_bbox_max = tiled_network.bbox;
        tiled_network.translate(bbox_min - tiled_bbox_min);
        return tiled_network;

    def __tile_with_mesh(self, mesh):
        tiler = WireTiler();
        tiler.set_single_cell_from_wire_network(self.unit_pattern);
        tiler.tile_hex_mesh(mesh, self.parameters);
        tiled_network = tiler.wire_network;
        return tiled_network;

    def __inflate(self, wire_network, options):
        if options.get("trim", False):
            wire_network.trim();

        inflator = PyWiresInflator(wire_network, self.parameters,
                options.get("periodic", False));
        inflator.inflate(
                subdivide_method="loop",
                subdivide_order = options.get("subdiv", 1));
        return inflator.mesh;

