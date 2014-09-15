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
        XY= np.array([-1,-1, 1]);
        YZ= np.array([ 1,-1,-1]);
        ZX= np.array([-1, 1,-1]);
        XYZ=np.array([-1,-1,-1]);

        self.symmetries = [
                lambda v: (v-bbox_center)*X + bbox_center,
                lambda v: (v-bbox_center)*Y + bbox_center,
                lambda v: (v-bbox_center)*Z + bbox_center,
                lambda v: (v-bbox_center)*XY + bbox_center,
                lambda v: (v-bbox_center)*YZ + bbox_center,
                lambda v: (v-bbox_center)*ZX + bbox_center,
                lambda v: (v-bbox_center)*XYZ + bbox_center,
                ];

    @property
    def value(self):
        raise NotImplementedError("This attribute is abstract");
