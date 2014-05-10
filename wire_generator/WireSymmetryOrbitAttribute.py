import numpy as np
from WireAttribute import WireAttribute

import PyMeshSetting
import PyMesh

class WireSymmetryOrbitAttribute(WireAttribute):
    def compute(self, wire_network):
        self.wire_network = wire_network;

        bbox_min, bbox_max = wire_network.bbox;
        bbox_center = 0.5 * (bbox_min + bbox_max);
        bbox_size = bbox_max - bbox_min;

        self.__initialize_grid(0.001 * np.amax(bbox_size));
        self.__initialize_reflective_symmetries(bbox_center);
        self.__compute_orbits();

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

    def __compute_orbits(self):
        self.orbits = np.arange(len(self.wire_network.vertices), dtype=int);
        for symm_map in self.symmetries:
            self.__collect_orbits(symm_map);

        indices = set(self.orbits);
        for i, key in enumerate(indices):
            self.orbits[self.orbits == key] = i;

    def __collect_orbits(self, symmetry_map):
        for i,v in enumerate(self.wire_network.vertices):
            v_reflected= symmetry_map(v);
            mapped_candidates = self.grid.get_items_near_point(v_reflected).ravel();
            if len(mapped_candidates) == 0:
                continue;
            vertex_group = [i] + mapped_candidates.tolist();
            self.orbits[vertex_group] = np.amin(self.orbits[vertex_group]);

    @property
    def value(self):
        return self.orbits

    @value.setter
    def value(self, val):
        raise RuntimeError("Wire symmetric orbit attribute is read only!");
