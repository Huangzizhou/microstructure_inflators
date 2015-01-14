import numpy as np
from InflatorFacade import InflatorFacade
from WireTiler import WireTiler
from WireInflator import WireInflator
from PeriodicWireInflator import PeriodicWireInflator
from PyWiresInflator import PyWiresInflator
from PyWiresTiler import PyWiresTiler

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
        repetitions = np.array(repetitions, dtype=int);

        tiler = PyWiresTiler();
        tiler.set_single_cell_from_wire_network(self.unit_pattern);
        tiler.tile(bbox_min, bbox_max, repetitions, self.parameters);

        tiled_network = tiler.wire_network;

        tiled_bbox_min, tiled_bbox_max = tiled_network.bbox;
        tiled_network.translate(bbox_min - tiled_bbox_min);
        return tiled_network;

    def __tile_with_mesh(self, mesh):
        tiler = PyWiresTiler();
        tiler.set_single_cell_from_wire_network(self.unit_pattern);
        tiler.tile_with_hex_mesh(mesh, self.parameters);
        tiled_network = tiler.wire_network;
        return tiled_network;

    def __inflate(self, wire_network, options):
        if options.get("trim", False):
            wire_network.trim();

        inflator = PyWiresInflator(wire_network, self.parameters,
                options.get("periodic", False));
        inflator.inflate(
                subdivide_method=options.get("subdiv_method", "simple"),
                subdivide_order = options.get("subdiv", 1),
                rel_geometry_correction = options.get("rel_geometry_correction"),
                abs_geometry_correction = options.get("abs_geometry_correction"),
                geometry_correction_cap = options.get("geometry_correction_cap"),
                geometry_spread = options.get("geometry_spread"),
                geometry_correction_lookup = options.get("geometry_correction_lookup"));
        return inflator.mesh;

