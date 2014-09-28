import numpy as np
from numpy.linalg import norm
import os
import os.path
import core.PyWireInflator2DSetting
from wire_io.WireWriter import WireWriter
import PyMesh
import PyWireInflator2D

from ParameterHandler2D import ParameterHandler2D

class WireInflator2D(object):
    def __init__(self, wire_network, parameters):
        self.wire_network = wire_network;
        self.__create_2D_inflator();
        self.parameters = parameters;
        self.parameter_handler = ParameterHandler2D(
                self.wire_network, self.inflator);

    def tile_periodic(self, bbox_min, bbox_max):
        self.__compute_scale_factor(bbox_min, bbox_max, 1, 1);
        self.inflator.set_dimension(1, 1);
        p = self.parameter_handler.convert_to_flattened_parameters(
                self.parameters);
        self.__scale_thickness_parameters(p);
        self.inflator.set_parameter(0, 0, p);
        self.inflator.set_max_triangle_area(0.0001);
        self.inflator.generate_periodic_pattern();
        self.is_periodic = True;

    def tile_box(self, bbox_min, bbox_max, rows, cols):
        self.__compute_scale_factor(bbox_min, bbox_max, rows, cols);
        self.inflator.set_dimension(rows, cols);
        for i in range(rows):
            for j in range(cols):
                p = self.parameter_handler.convert_to_flattened_parameters(
                        self.parameters);
                self.__scale_thickness_parameters(p);
                self.inflator.set_parameter(i, j, p);

        self.inflator.set_max_triangle_area(0.001);
        self.inflator.generate_tiled_pattern();

    def tile_quad_mesh(self, quad_mesh):
        raise NotImplementedError("Tiling with quad mesh is not supported.");

    def __create_2D_inflator(self):
        tmp_dir = "/tmp";
        path,name = os.path.split(self.wire_network.source_file);
        tmp_wire_file = os.path.join(tmp_dir, name);
        dim = self.wire_network.dim;

        if dim == 2:
            vertices = np.hstack((self.wire_network.vertices,
                np.zeros((self.wire_network.num_vertices, 1))));
            edges = self.wire_network.edges;

            writer = WireWriter(tmp_wire_file);
            writer.write(vertices, edges);

            self.inflator = PyWireInflator2D.WireInflatorFacade(tmp_wire_file);
            os.remove(tmp_wire_file);
        else:
            raise NotImplementedError("Wire network is not 2D.");

    def __compute_scale_factor(self, bbox_min, bbox_max, rows, cols):
        tol = 1e-3;
        dimension = np.array([cols, rows], dtype=float);
        scale_factor = np.divide(bbox_max - bbox_min, dimension);
        self.scale_factor = np.mean(scale_factor);
        if norm(scale_factor - self.scale_factor) > tol:
            raise RuntimeError("Non-uniform scaling {} is not supported!".format(
                scale_factor));

    def __scale_thickness_parameters(self, p):
        for i,value in enumerate(p):
            if self.inflator.get_parameter_type(i) ==\
                    PyWireInflator2D.WireInflatorFacade.THICKNESS:
                        # input: thickness in with unit mm.
                        # output: ratio of radius (thickness * 0.5) to cell
                        # size.  In other words, the wire radius if we scale the
                        # cell into a unit box.
                        p[i] = value / self.scale_factor * 0.5;

    @property
    def mesh(self):
        vertices = self.inflator.get_vertices().reshape((-1, 2), order="C");
        vertices = vertices * self.scale_factor;
        faces = self.inflator.get_triangles().reshape((-1, 3), order="C");

        factory = PyMesh.MeshFactory();
        factory.load_data(
                vertices.ravel(order="C"),
                faces.ravel(order="C"),
                np.zeros(0),
                2, 3, 0);
        mesh = factory.create_shared();

        attr_names = [];
        if hasattr(self, "is_periodic") and self.is_periodic:
            num_vertices = mesh.get_num_vertices();
            mesh.add_attribute("vertex_normal");
            normals = mesh.get_attribute("vertex_normal").reshape((num_vertices, -1));
            velocity = self.inflator.get_boundary_velocity();
            num_params = velocity.shape[1];
            for i in range(num_params):
                attr_name = "normal_velocity_{}".format(i);
                attr_value = normals * velocity[:,i][:,np.newaxis];
                mesh.add_attribute(attr_name);
                mesh.set_attribute(attr_name, attr_value.ravel(order="C"));
                attr_names.append(attr_name);
        return mesh;

