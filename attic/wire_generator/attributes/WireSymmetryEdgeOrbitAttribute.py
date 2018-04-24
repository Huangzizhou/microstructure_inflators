import numpy as np
from WireSymmetryAttribute import WireSymmetryAttribute

class WireSymmetryEdgeOrbitAttribute(WireSymmetryAttribute):
    def compute(self, wire_network):
        super(WireSymmetryEdgeOrbitAttribute, self).compute(wire_network);
        self.__generate_edge_lookup_table();
        self.__compute_vertex_symmetry_maps();
        self.__compute_edge_symmetries();

    def __generate_edge_lookup_table(self):
        self.edge_table = {};
        for i,edge in enumerate(self.wire_network.edges):
            e = tuple(sorted(edge));
            self.edge_table[e] = i;

    def __compute_vertex_symmetry_maps(self):
        self.vertex_maps = [];
        for symm_map in self.symmetries:
            v_map = self.__collect_orbits(symm_map);
            self.vertex_maps.append(v_map);

    def __collect_orbits(self, symmetry_map):
        v_map = np.zeros(self.wire_network.num_vertices, dtype=int) - 1;
        for i,v in enumerate(self.wire_network.vertices):
            if v_map[i] != -1:
                continue;

            v_reflected= symmetry_map(v);
            mapped_candidates = self.grid.get_items_near_point(v_reflected).ravel();
            if len(mapped_candidates) == 0:
                continue;
            elif len(mapped_candidates) != 1:
                raise RuntimeError("More than vertex mapped by symmetry!" +
                        "  Check for duplicated vertices.");

            j = mapped_candidates[0];
            v_map[i] = j;
            v_map[j] = i;
        return v_map;

    def __compute_edge_symmetries(self):
        edge_orbits = np.arange(self.wire_network.num_edges, dtype=int);
        for v_map in self.vertex_maps:
            self.__collect_edge_orbits(v_map, edge_orbits);

        edge_orbit_set = set(edge_orbits);
        for i,orbit_id in enumerate(edge_orbit_set):
            edge_orbits[edge_orbits == orbit_id] = i;

        self.edge_orbits = edge_orbits;

    def __collect_edge_orbits(self, v_map, edge_orbits):
        for i,ei in enumerate(self.wire_network.edges):
            ej = tuple(sorted(v_map[ei]));
            if ej not in self.edge_table:
                continue;
            j = self.edge_table[ej];

            orbit_i = edge_orbits[i];
            orbit_j = edge_orbits[j];

            self.__merge_orbits(edge_orbits, orbit_i, orbit_j);

    def __merge_orbits(self, orbits, orbit_i, orbit_j):
        idx = min(orbit_i, orbit_j);
        orbits[orbits == orbit_i] = idx;
        orbits[orbits == orbit_j] = idx;

    @property
    def value(self):
        return self.edge_orbits;
