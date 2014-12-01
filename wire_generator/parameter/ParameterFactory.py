import json
import numpy as np

from VertexThicknessParameter import VertexThicknessParameter
from EdgeThicknessParameter import EdgeThicknessParameter
from VertexOffsetParameter import VertexOffsetParameter

class ParameterFactory(object):
    def __init__(self, wire_network, default_thickness=0.5):
        self.wire_network = wire_network;
        self.default_thickness = default_thickness;

    def create_parameters_from_file(self, modifier_file=None):
        config = self.__load_modifier(modifier_file);
        self.create_parameters_from_dict(config);

    def create_parameters_from_dict(self, config):
        self.vertex_thickness_parameters = [];
        self.edge_thickness_parameters = [];
        self.vertex_offset_parameters = [];
        self.orbit_type = config.get("orbit_type", "orthotropic");

        if "thickness" in config:
            thickness_config = config["thickness"];
            thickness_type = thickness_config["type"];
            self.default_thickness = thickness_config["default"];
            if thickness_type == "vertex_orbit":
                self.__create_default_vertex_thickness_parameters();
            elif thickness_type == "edge_orbit":
                self.__create_default_edge_thickness_parameters();
            else:
                raise NotImplementedError(
                        "Unknown thickness type: {}".format(thickness_type));
            self.__modify_thickness_parameters(thickness_config);
        else:
            self.__create_default_vertex_thickness_parameters();

        if "vertex_offset" in config:
            offset_config = config["vertex_offset"];
            offset_type = offset_config["type"];
            assert(offset_type == "vertex_orbit");

            self.__create_default_vertex_offset_parameters();
            self.__modify_offset_parameters(offset_config);

    def __load_modifier(self, modifier_file):
        if modifier_file is None: return {}
        with open(modifier_file, 'r') as fin:
            config = json.load(fin);
        return config;

    def __create_default_vertex_thickness_parameters(self):
        num_orbits = self.__get_num_vertex_orbits();
        self.vertex_thickness_parameters = [
                VertexThicknessParameter(self.wire_network, i, self.default_thickness)
                for i in range(num_orbits) ];

        for param in self.vertex_thickness_parameters:
            param.orbit_type = self.orbit_type;

    def __create_default_edge_thickness_parameters(self):
        num_orbits = self.__get_num_edge_orbits();
        self.edge_thickness_parameters = [
                EdgeThicknessParameter(self.wire_network, i, self.default_thickness)
                for i in range(num_orbits) ];

        for param in self.edge_thickness_parameters:
            param.orbit_type = self.orbit_type;

    def __create_default_vertex_offset_parameters(self):
        num_orbits = self.__get_num_vertex_orbits();
        self.vertex_offset_parameters = [
                VertexOffsetParameter(self.wire_network, i)
                for i in range(num_orbits) ];

        for param in self.vertex_offset_parameters:
            param.orbit_type = self.orbit_type;

    def __modify_thickness_parameters(self, thickness_config):
        effective_orbits = thickness_config["effective_orbits"];
        effective_values = thickness_config["thickness"];
        if thickness_config["type"] == "vertex_orbit":
            thickness_params = self.vertex_thickness_parameters;
        else:
            thickness_params = self.edge_thickness_parameters;

        for orbit_id, thickness in zip(effective_orbits, effective_values):
            param = thickness_params[orbit_id];
            assert(param.orbit_id == orbit_id);
            param.set_formula(thickness);

    def __modify_offset_parameters(self, offset_config):
        effective_orbits = offset_config["effective_orbits"];
        effective_values = offset_config["offset_percentages"];
        for orbit_id, offset in zip(effective_orbits, effective_values):
            param = self.vertex_offset_parameters[orbit_id];
            assert(param.orbit_id == orbit_id);
            param.set_formula(offset);

    def __get_num_vertex_orbits(self):
        if self.orbit_type == "isotropic":
            return self.__get_num_isotropic_vertex_orbits();
        elif self.orbit_type == "orthotropic":
            return self.__get_num_orthotropic_vertex_orbits();
        else:
            raise NotImplementedError("Unsupported orbit type: {}"\
                    .format(self.orbit_type));

    def __get_num_edge_orbits(self):
        if self.orbit_type == "isotropic":
            return self.__get_num_isotropic_edge_orbits();
        elif self.orbit_type == "orthotropic":
            return self.__get_num_orthotropic_edge_orbits();
        else:
            raise NotImplementedError("Unsupported orbit type: {}"\
                    .format(self.orbit_type));

    def __get_num_orthotropic_vertex_orbits(self):
        if "orthotropic_symmetry_vertex_orbit" not in self.wire_network.attributes:
            self.wire_network.compute_symmetry_orbits();
        return len(np.unique(
            self.wire_network.attributes["orthotropic_symmetry_vertex_orbit"]));

    def __get_num_isotropic_vertex_orbits(self):
        if "isotropic_symmetry_vertex_orbit" not in self.wire_network.attributes:
            self.wire_network.compute_symmetry_orbits();
        return len(np.unique(
            self.wire_network.attributes["isotropic_symmetry_vertex_orbit"]));

    def __get_num_isotropic_edge_orbits(self):
        if "isotropic_symmetry_edge_orbit" not in self.wire_network.attributes:
            self.wire_network.compute_symmetry_orbits();
        return len(np.unique(self.wire_network.attributes["isotropic_symmetry_edge_orbit"]));

    def __get_num_orthotropic_edge_orbits(self):
        if "orthotropic_symmetry_edge_orbit" not in self.wire_network.attributes:
            self.wire_network.compute_symmetry_orbits();
        return len(np.unique(self.wire_network.attributes["orthotropic_symmetry_edge_orbit"]));

    def is_initialized(self):
        return hasattr(self, "vertex_thickness_parameters") and\
                hasattr(self, "edge_thickness_parameters") and\
                hasattr(self, "vertex_offset_parameters")

    @property
    def parameters(self):
        if not self.is_initialized():
            raise RuntimeError("Parameters are not initialized.");

        return self.vertex_thickness_parameters + \
                self.edge_thickness_parameters + \
                self.vertex_offset_parameters;

