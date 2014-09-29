import numpy as np
from InflatorFacade import InflatorFacade
from WireInflator2D import WireInflator2D

class InflatorFacade2D(InflatorFacade):
    def __init__(self, wire_network, parameters):
        super(InflatorFacade2D, self).__init__(wire_network, parameters);

    def inflate_with_guide_box(self, bbox_min, bbox_max, repetitions, options):
        inflator = WireInflator2D(self.unit_pattern, self.parameters);
        if options.get("trim", False):
            raise NotImplementedError("Trimming is not supported in 2D");

        bbox_min = np.array(bbox_min)[:2];
        bbox_max = np.array(bbox_max)[:2];
        if options.get("periodic", False):
            inflator.tile_periodic(bbox_min, bbox_max);
        else:
            cols, rows = repetitions[:2];
            inflator.tile_box(bbox_min, bbox_max, rows, cols);

        return inflator.mesh;

    def inflate_with_guide_mesh(self, mesh, options):
        inflator = WireInflator2D(self.unit_pattern, self.parameters);
        if options.get("trim", False):
            raise NotImplementedError("Trimming is not supported in 2D");
        if options.get("periodic", False):
            raise RuntimeError(
                    "Perodic inflation is not valid when tiling with guide mesh");

        inflator.tile_quad_mesh(mesh);
        return inflator.mesh;
