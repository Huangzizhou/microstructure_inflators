import numpy as np
from WireAttribute import WireAttribute

import PyMeshSetting
import PyMesh

class WireSymmetryAttribute(WireAttribute):
    def compute(self, wire_network):
        self.wire_network = wire_network;

        bbox_min, bbox_max = wire_network.bbox;
        bbox_center = 0.5 * (bbox_min + bbox_max);
        bbox_size = bbox_max - bbox_min;

        self.__initialize_grid(0.001 * np.amax(bbox_size));
        self.__initialize_reflective_symmetries(bbox_center);

    def __initialize_grid(self, cell_size):
        self.grid = PyMesh.HashGrid.create(cell_size);
        for i,v in enumerate(self.wire_network.vertices):
            self.grid.insert(i, v);

    def __initialize_reflective_symmetries(self, bbox_center):
        X = np.array([-1, 1, 1]);
        Y = np.array([ 1,-1, 1]);
        Z = np.array([ 1, 1,-1]);
        reflect_X = lambda v: (v-bbox_center) * X + bbox_center
        reflect_Y = lambda v: (v-bbox_center) * Y + bbox_center
        reflect_Z = lambda v: (v-bbox_center) * Z + bbox_center
        self.symmetries = [reflect_X, reflect_Y, reflect_Z];

    @property
    def value(self):
        raise NotImplementedError("This attribute is abstract");
