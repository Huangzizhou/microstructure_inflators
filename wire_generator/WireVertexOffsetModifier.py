from WireModifier import WireModifier

import json
import numpy as np

class WireVertexOffsetModifier(WireModifier):
    def __init__(self, config):
        """ Syntax:
        "vertex_offset": {
            "type": "vertex_orbit",
            "effective_orbits": [i0, i1, ...],
            "offset_percentages": [[#, #, #], [0.0, #, 0.0], ...,\
                    ["{x}", "{y}", 0.0]]
        }
        """
        self.config = config;
        self.offset_type = config["type"];
        self.offset_percentages = config["offset_percentages"];

    def modify(self, wire_network, **kwargs):
        self.__load_orbits(wire_network);
        offsets = self.__generate_default_offsets(wire_network);
        self.__compute_offset(wire_network, offsets, **kwargs);
        self.__assign_offset_attribute(wire_network, offsets);

    def __load_orbits(self, wire_network):
        if "symmetry_vertex_orbit" not in wire_network.attributes:
            wire_network.compute_symmetry_orbits();

        if self.offset_type == "vertex_orbit":
            self.orbits = wire_network.attributes["symmetry_vertex_orbit"];
        else:
            raise NotImplementedError("Vertex offset type ({}) is not supported"\
                    .format(self.offset_type));

        self.orbits = np.array(self.orbits);
        self.effective_orbits = self.orbits[self.config["effective_orbits"]];

    def __generate_default_offsets(self, wire_network):
        offsets = np.zeros((wire_network.num_vertices, wire_network.dim));
        return offsets;

    def __compute_offset(self, wire_network, offsets, **kwargs):
        for orbit, offset_percent in zip(
                self.effective_orbits, self.offset_percentages):
            offset = self.__compute_offset_per_orbit(wire_network, orbit,
                    offset_percent, **kwargs);
            offsets[orbit] = offset;

    def __compute_offset_per_orbit(
            self, wire_network, orbit, offset_percent, **kwargs):
        vertices = wire_network.vertices[orbit];
        bbox_min = np.amin(vertices, axis=0);
        bbox_max = np.amax(vertices, axis=0);
        bbox_center = 0.5 * (bbox_min + bbox_max);
        bbox_size = bbox_max - bbox_min;

        offset_percent = [eval(p.format(**kwargs)) if isinstance(p, str) else p
                for p in offset_percent];

        offset = [(v-bbox_center)*offset_percent for v in vertices];
        return offset;

    def __assign_offset_attribute(self, wire_network, offsets):
        attr_name = "vertex_offset";
        if attr_name not in wire_network.attributes:
            wire_network.attributes.add(attr_name, offsets);
        else:
            wire_network.attributes[attr_name] = offsets;

