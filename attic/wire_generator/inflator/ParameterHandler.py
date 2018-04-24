import numpy as np
import re
from parameter.VertexThicknessParameter import VertexThicknessParameter
from parameter.EdgeThicknessParameter import EdgeThicknessParameter
from parameter.VertexOffsetParameter import VertexOffsetParameter

import PyWires

class ParameterHandler(object):
    def __init__(self, wire_network):
        raise DeprecationWarning("This method is deprecated.");
        self.wire_network = wire_network;
        self.__remove_old_attributes();
        self.__initialize_orbits();

    def __remove_old_attributes(self):
        if "vertex_thickness" in self.wire_network.attributes:
            del self.wire_network.attributes["vertex_thickness"];
        if "edge_thickness" in self.wire_network.attributes:
            del self.wire_network.attributes["edge_thickness"];
        if "vertex_offset" in self.wire_network.attributes:
            del self.wire_network.attributes["vertex_offset"];

    def __initialize_orbits(self):
        if "orthotropic_symmetry_vertex_orbit" not in self.wire_network.attributes or\
                "orthotropic_symmetry_edge_orbit" not in self.wire_network.attributes or\
                "isotropic_symmetry_edge_orbit" not in self.wire_network.attributes or\
                "isotropic_symmetry_vertex_orbit" not in self.wire_network.attributes:
            self.wire_network.compute_symmetry_orbits();

        self.orthotropic_vertex_orbits = \
                self.wire_network.attributes["orthotropic_symmetry_vertex_orbit"];
        self.isotropic_vertex_orbits = \
                self.wire_network.attributes["isotropic_symmetry_vertex_orbit"];
        self.orthotropic_edge_orbits = \
                self.wire_network.attributes["orthotropic_symmetry_edge_orbit"];
        self.isotropic_edge_orbits = \
                self.wire_network.attributes["isotropic_symmetry_edge_orbit"];

    def convert_to_attributes(self, parameters, **kwargs):
        attributes = self.wire_network.attributes;
        for param in parameters:
            if param.orbit_type == "isotropic":
                vertex_orbits = self.isotropic_vertex_orbits;
                edge_orbits = self.isotropic_edge_orbits;
            elif param.orbit_type == "orthotropic":
                vertex_orbits = self.orthotropic_vertex_orbits;
                edge_orbits = self.isotropic_edge_orbits;
            else:
                raise NotImplementedError("Orbit type \"{}\" is not supported"\
                        .format(param.orbit_type));

            if isinstance(param, VertexThicknessParameter):
                if "vertex_thickness" not in attributes:
                    attributes["vertex_thickness"] = np.zeros(
                            self.wire_network.num_vertices);
                thickness = attributes["vertex_thickness"];
                affected_vertices = vertex_orbits == param.orbit_id;
                thickness[affected_vertices] = param.evaluate(**kwargs);
            elif isinstance(param, EdgeThicknessParameter):
                if "edge_thickness" not in attributes:
                    attributes["edge_thickness"] = np.zeros(
                            self.wire_network.num_edges);
                thickness = attributes["edge_thickness"];
                affected_edges = edge_orbits == param.orbit_id;
                thickness[affected_edges] = param.evaluate(**kwargs);
            elif isinstance(param, VertexOffsetParameter):
                if "vertex_offset" not in attributes:
                    attributes["vertex_offset"] = np.zeros((
                        self.wire_network.num_vertices,
                        self.wire_network.dim));
                offset = attributes["vertex_offset"];
                affected_vertices = vertex_orbits == param.orbit_id;
                percentage = param.evaluate(**kwargs);
                vertices = self.wire_network.vertices[affected_vertices];
                centroid = self.wire_network.bbox_center;
                offset[affected_vertices] = (vertices - centroid) * percentage;
            else:
                raise NotImplementedError(
                        "Unknown parameter type: {}".format(param));

    def convert_to_PyWires_parameter_manager(self, parameters):
        default_thickness = 0.5;
        parameter_manager = PyWires.ParameterManager.create_empty_manager(
                self.wire_network.raw_wires, default_thickness);
        thickness_type = None;

        vertex_indices = np.arange(self.wire_network.num_vertices, dtype=int);
        edge_indices = np.arange(self.wire_network.num_edges, dtype=int);

        for param in parameters:
            if param.orbit_type == "isotropic":
                vertex_orbits = self.isotropic_vertex_orbits;
                edge_orbits = self.isotropic_edge_orbits;
            elif param.orbit_type == "orthotropic":
                vertex_orbits = self.orthotropic_vertex_orbits;
                edge_orbits = self.isotropic_edge_orbits;
            else:
                raise NotImplementedError("Orbit type \"{}\" is not supported"\
                        .format(param.orbit_type));

            if isinstance(param, VertexThicknessParameter):
                if thickness_type is None:
                    thickness_type = PyWires.VERTEX;
                    parameter_manager.set_thickness_type(thickness_type);
                else:
                    assert(thickness_type == PyWires.VERTEX);

                roi = vertex_indices[vertex_orbits == param.orbit_id];
                formula, value = self.__separate_formula_and_value(param,
                        param.default_thickness);
                parameter_manager.add_thickness_parameter(roi, formula, value);
            elif isinstance(param, EdgeThicknessParameter):
                if thickness_type is None:
                    thickness_type = PyWires.EDGE;
                    parameter_manager.set_thickness_type(thickness_type);
                else:
                    assert(thickness_type == PyWires.EDGE);

                roi = edge_indices[edge_orbits == param.orbit_id];
                formula, value = self.__separate_formula_and_value(param,
                        param.default_thickness);
                parameter_manager.add_thickness_parameter(roi, formula, value);
            elif isinstance(param, VertexOffsetParameter):
                roi = vertex_indices[vertex_orbits == param.orbit_id];
                offset_percentage = param.default_offset;
                if hasattr(param, "formula"):
                    offset_percentage = param.formula;

                for i,entry in enumerate(offset_percentage):
                    formula, value = self.__split_formula(entry);
                    parameter_manager.add_offset_parameter(
                            roi, formula, value, i);
            else:
                raise NotImplementedError(
                        "Unknown parameter type: {}".format(param));

        return parameter_manager;

    def __separate_formula_and_value(self, param, default_value=0.0):
        if hasattr(param, "formula"):
            return self.__split_formula(param.formula, default_value);
        else:
            return "", default_value;

    def __split_formula(self, data, default_value=0.0):
        if isinstance(data, (str, unicode)):
            formula = str(data);
            result = re.match("{(.*)}", formula);
            if (result is not None):
                formula = result.group(1);
            else:
                raise RuntimeError("Invalid formula: \"{}\"".format(formula));
            value = default_value;
        elif isinstance(data, float):
            formula = "";
            value = data;
        return formula, value;

