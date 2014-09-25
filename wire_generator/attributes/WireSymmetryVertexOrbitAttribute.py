import numpy as np
from WireSymmetryAttribute import WireSymmetryAttribute

import PyMesh

class WireSymmetryVertexOrbitAttribute(WireSymmetryAttribute):
    def compute(self, wire_network):
        super(WireSymmetryVertexOrbitAttribute, self).compute(wire_network);
        self.__compute_orbits();

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
