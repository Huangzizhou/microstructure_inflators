from WireModifier import WireModifier

import json
import numpy as np

class WireVertexOffsetModifier(WireModifier):
    def __init__(self, config):
        """ Syntax:
        {
            "type": "vertex_orbit",
            "orbit_file": "filename",
            "effective_orbits": [i0, i1, ...],
            "offset_percentages": [[#, #, #], [0.0, #, 0.0], ...,\
                    ["{x}", "{y}", 0.0]]
        }
        """
        self.offset_type = config["type"];
        self.__load_orbits(config["orbit_file"]);
        self.effective_orbits = self.orbits[config["effective_orbits"]];
        self.offset_percentages = config["offset_percentages"];

    def modify(self, wire_network, **kwargs):
        offsets = self.__generate_default_offsets(wire_network);
        self.__compute_offset(wire_network, offsets, **kwargs);
        #self.__update_wire_network(wire_network, offsets);
        self.__assign_offset_attribute(wire_network, offsets);

    def __load_orbits(self, orbit_file):
        with open(orbit_file, 'r') as fin:
            contents = json.load(fin);
            if self.offset_type == "vertex_orbit":
                orbits = contents["vertex_orbits"];
            else:
                raise NotImplementedError("Vertex offset type ({}) is not supported"\
                        .format(self.offset_type));
            self.orbits = np.array(orbits);

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

    def __update_wire_network(self, wire_network, offsets):
        for i in range(wire_network.num_vertices):
            wire_network.vertices[i] += offsets[i];

    def __assign_offset_attribute(self, wire_network, offsets):
        attr_name = "vertex_offset";
        if attr_name not in wire_network.attributes:
            wire_network.attributes.add(attr_name, offsets);
        else:
            wire_network.attributes[attr_name] = offsets;

