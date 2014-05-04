import unittest
import numpy as np

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from MaterialFitterFactory import MaterialFitterFactory

class MaterialFitterTest(unittest.TestCase):
    def add_attribute_to_mesh(self, mesh, attr_name, attr_value):
        if mesh.has_attribute(attr_name):
            raise RuntimeError("Attribute {} already exists.".format(attr_name));
        mesh.add_attribute(attr_name);
        mesh.set_attribute(attr_name, attr_value);

    def compute_vertex_voronoi(self, mesh):
        if not mesh.has_attribute("face_voronoi_area"):
            mesh.add_attribute("face_voronoi_area");
        face_voronoi_areas = mesh.get_face_attribute("face_voronoi_area");
        vertex_voronoi_areas = np.zeros(mesh.num_vertices);
        on_face = np.zeros(mesh.num_vertices, dtype=bool);
        for i,area in enumerate(face_voronoi_areas):
            face = mesh.faces[i];
            on_face[face[0]] = True;
            on_face[face[1]] = True;
            on_face[face[2]] = True;
            vertex_voronoi_areas[face] += area;

        mesh.add_attribute("vertex_voronoi_area");
        mesh.set_attribute("vertex_voronoi_area", vertex_voronoi_areas);

    def save_mesh_fields(self, mesh_file, mesh, pressures=None, displacements=None,
            stress_traces=None):
        num_fields = len(displacements);
        attributes_to_save = [];
        for i in range(num_fields):
            if pressures is not None:
                pressure_attr_name = "pressure_{}".format(i);
                self.add_attribute_to_mesh(mesh, pressure_attr_name, pressures[i]);
                attributes_to_save.append(pressure_attr_name);

            if displacements is not None:
                disp_attr_name = "displacement_{}".format(i);
                self.add_attribute_to_mesh(mesh, disp_attr_name, displacements[i]);
                attributes_to_save.append(disp_attr_name);

            if stress_traces is not None:
                strs_attr_name = "stress_trace_{}".format(i);
                self.add_attribute_to_mesh(mesh, strs_attr_name, stress_traces[i]);
                attributes_to_save.append(strs_attr_name);
        self.compute_vertex_voronoi(mesh);
        attributes_to_save.append("vertex_voronoi_area");
        save_mesh(mesh_file, mesh, *attributes_to_save);

    def assertEqualMaterial(self, dim, mat1, mat2):
        ori = np.zeros(dim);
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for l in range(dim):
                        self.assertAlmostEqual(
                                mat1.get_material_tensor(i,j,k,l,ori),
                                mat2.get_material_tensor(i,j,k,l,ori));

