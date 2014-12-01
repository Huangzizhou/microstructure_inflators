import numpy as np

from parameter.VertexThicknessParameter import VertexThicknessParameter
from parameter.EdgeThicknessParameter import EdgeThicknessParameter
from parameter.VertexOffsetParameter import VertexOffsetParameter

import core.PyWireInflator2DSetting
import PyWireInflator2D

class ParameterHandler2D(object):
    def __init__(self, wire_network, inflator):
        self.wire_network = wire_network;
        self.inflator = inflator;
        self.__initialize_orbits();
        self.__compute_orbit_index_maps();

    def __initialize_orbits(self):
        if "orthotropic_symmetry_vertex_orbit" not in self.wire_network.attributes or\
                "isotropic_symmetry_vertex_orbit" not in self.wire_network.attributes or\
                "orthotropic_symmetry_edge_orbit" not in self.wire_network.attributes or\
                "isotropic_symmetry_edge_orbit" not in self.wire_network.attributes:
            self.wire_network.compute_symmetry_orbits();

        self.orthotropic_vertex_orbits = \
                self.wire_network.attributes["orthotropic_symmetry_vertex_orbit"];
        self.isotropic_vertex_orbits = \
                self.wire_network.attributes["isotropic_symmetry_vertex_orbit"];
        self.orthotropic_edge_orbits = \
                self.wire_network.attributes["orthotropic_symmetry_edge_orbit"];
        self.isotropic_edge_orbits = \
                self.wire_network.attributes["isotropic_symmetry_edge_orbit"];

    def __compute_orbit_index_maps(self):
        num_parameters = self.inflator.get_num_parameters();

        # orbit maps:
        # thickness: {vertex indices} -> parameter index
        # offset:    {vertex indices} -> [parameter index, parameter index]
        self.thickness_orbit_map = {};
        self.offset_orbit_map = {};

        for i in range(num_parameters):
            param_type = self.inflator.get_parameter_type(i);
            if param_type == PyWireInflator2D.WireInflatorFacade.VERTEX_OFFSET:
                self.__register_offset_orbit(i);
            elif param_type == PyWireInflator2D.WireInflatorFacade.THICKNESS:
                self.__register_thickness_orbit(i);
            else:
                raise NotImplementedError(
                        "Unknown parameter type: {}".format(param_type));

    def __register_offset_orbit(self, index):
        tol = 1e-3;
        affected_vertices = self.inflator.get_affected_vertex_orbit(index).ravel();
        affectec_vertices = sorted(affected_vertices.ravel().astype(int));
        offset_dir = self.inflator.get_offset_direction(index);
        mask = np.amax(np.absolute(offset_dir), axis=0) > tol;

        key = tuple(affected_vertices);
        dof_index = self.offset_orbit_map.get(key, np.ones(2, dtype=int) * -1);
        dof_index[mask] = index;
        self.offset_orbit_map[key] = dof_index;

    def __register_thickness_orbit(self, index):
        affected_vertices = self.inflator.get_affected_vertex_orbit(index).ravel();
        affectec_vertices = sorted(affected_vertices.ravel().astype(int));
        self.thickness_orbit_map[tuple(affected_vertices)] = index;

    def convert_to_flattened_parameters(self, parameters, **kwargs):
        vertex_indices = np.arange(self.wire_network.num_vertices, dtype=int);
        parameter = np.zeros(self.inflator.get_num_parameters());
        for param in parameters:
            if isinstance(param, VertexThicknessParameter):
                index = self.__lookup_thickness_vertex_orbit_map(param.orbit_id,
                        param.orbit_type);
                parameter[index] = param.evaluate(**kwargs);
            elif isinstance(param, EdgeThicknessParameter):
                raise NotImplementedError("Edge thickness is not supported in 2D");
            elif isinstance(param, VertexOffsetParameter):
                index = self.__lookup_offset_vertex_orbit_map(param.orbit_id,
                        param.orbit_type);
                mask = index >= 0;
                parameter[index[mask]] = param.evaluate(**kwargs)[mask];
            else:
                raise NotImplementedError(
                        "Unknown parameter type: {}".format(param));
        return parameter;

    def __lookup_thickness_vertex_orbit_map(self, orbit_id, orbit_type):
        vertex_indices = np.arange(self.wire_network.num_vertices, dtype=int);
        affected_vertices = self.__get_affected_vertices(orbit_id, orbit_type);
        key = tuple(vertex_indices[affected_vertices]);
        if key not in self.thickness_orbit_map:
            raise RuntimeError(
                    "Thickness parameter of vertex orbit {} is not valid.".format(orbit_id));
        return self.thickness_orbit_map[key];

    def __lookup_offset_vertex_orbit_map(self, orbit_id, orbit_type):
        vertex_indices = np.arange(self.wire_network.num_vertices, dtype=int);
        affected_vertices = self.__get_affected_vertices(orbit_id, orbit_type);
        key = tuple(vertex_indices[affected_vertices]);
        return self.offset_orbit_map.get(key, np.ones(2, dtype=int) * -1);

    def __get_affected_vertices(self, orbit_id, orbit_type):
        if orbit_type == "orthotropic":
            affected_vertices = self.orthotropic_vertex_orbits == orbit_id;
        elif orbit_type == "isotropic":
            affected_vertices = self.isotropic_vertex_orbits == orbit_id;
        else:
            raise RuntimeError("Unsupported orbit type: {}".format(orbit_type));
        return affected_vertices;

