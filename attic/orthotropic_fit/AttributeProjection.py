import numpy as np
from numpy.linalg import lstsq

import LinearElasticity
import PyMeshUtils
from Mesh import Mesh
from timethis import timethis

class AttributeProjection(object):
    @timethis
    def __init__(self, target_mesh):
        self.target_mesh = target_mesh;
        self.point_locator = PyMeshUtils.PointLocator(self.target_mesh.raw_mesh);

    @timethis
    def project(self, source_mesh, vertex_fields):
        assert(self.target_mesh.dim == source_mesh.dim);
        self.point_locator.locate(source_mesh.vertices);

        element_indices = self.point_locator.get_enclosing_voxels().ravel();
        barycentric_coords = self.point_locator.get_barycentric_coords();

        barycentric_matrix = np.zeros((source_mesh.num_vertices,
            self.target_mesh.num_vertices));
        for i in range(source_mesh.num_vertices):
            enclosing_elem = self.target_mesh.elements[element_indices[i]];
            barycentric_matrix[i, enclosing_elem] = barycentric_coords[i];

        projected_fields = [];
        for field in vertex_fields:
            field = field.reshape((source_mesh.num_vertices, -1), order="C");
            proj_field, residual, rank, singular_vals =\
                    lstsq(barycentric_matrix, field);
            projected_fields.append(proj_field.ravel(order="C"));

        return projected_fields;

