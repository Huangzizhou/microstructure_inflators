import json
import numpy as np
from numpy.linalg import norm

from WireModifier import WireModifier

import PyWireInflator2DSetting
import PyWireInflator2D

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
        assert(self.offset_type == "vertex_orbit");
        self.__load_orbits(config["orbit_file"]);
        self.effective_orbits = self.orbits[config["effective_orbits"]];
        self.offset_percentages = config["offset_percentages"];

    def modify(self, param, inflator, **kwargs):
        assert(len(param) == inflator.get_num_parameters());
        tol = 1e-3;
        self.__compute_offset_orbit_index_map(inflator);
        num_params = len(param);
        effective_offsets = self.__compute_offset(**kwargs);
        for orbit, offset in zip(self.effective_orbits, effective_offsets):
            orbit = tuple(sorted(orbit));
            param_index = self.orbit_index_map.get(orbit, None);
            if param_index is None:
                raise RuntimeError("Vertex orbit {} is not valid".format(
                    orbit));

            mask = np.absolute(offset) > tol
            if np.any(mask):
                param_index = param_index[mask];
                assert(param_index >= 0);
                param[param_index] = offset[mask];

        return param;

    def __load_orbits(self, orbit_file):
        with open(orbit_file, 'r') as fin:
            contents = json.load(fin);
            if self.offset_type == "vertex_orbit":
                orbits = contents["vertex_orbits"];
            elif self.offset_type == "edge_orbit":
                orbits = contents["edge_orbits"];
            else:
                raise NotImplementedError("Offset type ({}) is not supported"\
                        .format(self.offset_type));
            self.orbits = np.array(orbits);

    def __compute_offset_orbit_index_map(self, inflator):
        tol = 1e-3;
        target_param_type = PyWireInflator2D.WireInflatorFacade.VERTEX_OFFSET;
        num_parameters = inflator.get_num_parameters();
        self.orbit_index_map = {};
        for i in range(num_parameters):
            if inflator.get_parameter_type(i) != target_param_type:
                continue;
            vertex_orbit = sorted(inflator.get_affected_vertex_orbit(i).ravel().astype(int));
            offset_dir = inflator.get_offset_direction(i);
            offset_dir = np.amax(np.absolute(offset_dir), axis=0);
            mask = offset_dir > tol;
            key = tuple(vertex_orbit);
            val = self.orbit_index_map.get(key, np.ones(2, dtype=int)*-1);
            val[mask] = i;
            self.orbit_index_map[key] = val;

    def __compute_offset(self, **kwargs):
        effective_offsets = [];
        for orbit, offset_percent in zip(
                self.effective_orbits, self.offset_percentages):
            offset = self.__compute_offset_per_orbit(orbit, offset_percent, **kwargs);
            effective_offsets.append(offset);
        return np.array(effective_offsets);

    def __compute_offset_per_orbit(
            self, orbit, offset_percent, **kwargs):
        offset_percent = [eval(p.format(**kwargs))
                if isinstance(p, (str, unicode)) else p
                for p in offset_percent];
        return offset_percent;
        # TODO:
        # The generic vertex offset supports asymmetric offset in x and y
        # direction.  It essentially allows changing a square into a rectangle.
        # Currently, the 2D wire inflator only support symmetric offset, i.e.
        # only uniform scale of squares.  So this modifier convert asymmetric
        # offset into symmetric offset by averaging.
        #return np.mean(offset_percent);
