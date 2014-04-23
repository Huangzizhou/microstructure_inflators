import unittest
import numpy as np

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from OrthotropicFitter import OrthotropicFitter

def add_attribute_to_mesh(mesh, attr_name, attr_value):
    if mesh.has_attribute(attr_name):
        raise RuntimeError("Attribute {} already exists.".format(attr_name));
    mesh.add_attribute(attr_name);
    mesh.set_attribute(attr_name, attr_value);

def compute_vertex_voronoi(mesh):
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

def save_mesh_fields(mesh_file, mesh, pressures=None, displacements=None,
        stress_traces=None):
    num_fields = len(displacements);
    attributes_to_save = [];
    for i in range(num_fields):
        if pressures is not None:
            pressure_attr_name = "pressure_{}".format(i);
            add_attribute_to_mesh(mesh, pressure_attr_name, pressures[i]);
            attributes_to_save.append(pressure_attr_name);

        if displacements is not None:
            disp_attr_name = "displacement_{}".format(i);
            add_attribute_to_mesh(mesh, disp_attr_name, displacements[i]);
            attributes_to_save.append(disp_attr_name);

        if stress_traces is not None:
            strs_attr_name = "stress_trace_{}".format(i);
            add_attribute_to_mesh(mesh, strs_attr_name, stress_traces[i]);
            attributes_to_save.append(strs_attr_name);
    compute_vertex_voronoi(mesh);
    attributes_to_save.append("vertex_voronoi_area");
    save_mesh(mesh_file, mesh, *attributes_to_save);

class OrthotropicFitterTest(unittest.TestCase):
    def setUp(self):
        self.square = load_mesh("examples/solid_square.obj");
        self.cube = load_mesh("examples/solid_cube.msh");

    def test_2D(self):
        fitter = OrthotropicFitter(self.square, "examples/isotropic.material");
        fitter.fit();
        save_mesh_fields("tmp2D.msh", fitter.mesh,
                fitter.pressures,
                fitter.displacements,
                fitter.stress_traces);
        save_mesh_fields("tmp2D_coarse.msh", fitter.coarse_mesh,
                displacements = fitter.coarse_displacements);
        self.assertEqual(2, len(fitter.youngs_modulus));
        self.assertEqual(2, len(fitter.poisson_ratio));
        self.assertEqual(1, len(fitter.shear_modulus));
        self.assertAlmostEqual(1.0, fitter.youngs_modulus[0], delta=1e-3);
        self.assertAlmostEqual(1.0, fitter.youngs_modulus[1], delta=1e-3);
        self.assertAlmostEqual(0.0, fitter.poisson_ratio[0], delta=1e-3);
        self.assertAlmostEqual(0.0, fitter.poisson_ratio[1], delta=1e-3);
        self.assertAlmostEqual(0.5, fitter.shear_modulus[0], delta=1e-3);

    def test_3D(self):
        fitter = OrthotropicFitter(self.cube);
        fitter.fit();
        save_mesh_fields("tmp3D.msh", fitter.mesh,
                fitter.pressures,
                fitter.displacements,
                fitter.stress_traces);
        save_mesh_fields("tmp3D_coarse.msh", fitter.coarse_mesh,
                displacements = fitter.coarse_displacements);
        self.assertEqual(3, len(fitter.youngs_modulus));
        self.assertEqual(6, len(fitter.poisson_ratio));
        self.assertEqual(3, len(fitter.shear_modulus));
        self.assertAlmostEqual(1.0, np.amax(fitter.youngs_modulus), delta=1e-2);
        self.assertAlmostEqual(1.0, np.amin(fitter.youngs_modulus), delta=1e-2);
        self.assertAlmostEqual(0.0, np.amax(fitter.poisson_ratio), delta=1e-2);
        self.assertAlmostEqual(0.0, np.amin(fitter.poisson_ratio), delta=1e-2);
        self.assertAlmostEqual(0.5, np.amax(fitter.shear_modulus), delta=1e-2);
        self.assertAlmostEqual(0.5, np.amin(fitter.shear_modulus), delta=1e-2);

