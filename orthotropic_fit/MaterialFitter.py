import numpy as np
from numpy.linalg import lstsq, inv, norm

import LinearElasticitySettings
import PyMeshUtils

from Mesh import Mesh
from mesh_io import load_mesh, form_mesh
import PyAssembler

from AttributeProjection import AttributeProjection
from BoundaryConditionExtractor import BoundaryConditionExtractor
import ElasticModel2
from BoxMeshGenerator import generate_box_mesh
from ElasticityUtils import displacement_to_stress, displacement_to_strain,\
        total_energy
from LinearElasticity import LinearElasticity
from Material import Material
from timethis import timethis

class MaterialFitter(object):
    def __init__(self, mesh, material_file=None):
        self.mesh = mesh;
        mat = Material(self.mesh.dim, material_file);
        assembler = PyAssembler.FEAssembler.create(mesh.raw_mesh, mat.material);
        self.assembler = ElasticModel2.PyAssembler(mesh, assembler, mat);
        self.deformer = LinearElasticity(self.mesh, self.assembler);
        self.displacements = [];
        self.stress_traces = [];
        self.forces = [];
        self.pressures = [];

    @timethis
    def fit(self):
        self._generate_boundary_conditions();
        self._deform_shape();
        self._initialize_coarse_mesh();
        self._fit_material_parameters();

    @timethis
    def _generate_boundary_conditions(self):
        bd = BoundaryConditionExtractor(self.mesh);
        self.boundary_conditions = [];
        for config in self.bc_configs:
            bd.clear();
            bd.extract_from_dict(config);
            self.boundary_conditions.append([
                bd.neumann_bc, bd.dirichlet_bc]);

    @timethis
    def _deform_shape(self):
        dim = self.mesh.dim;
        stress_size = dim*(dim+1)/2;
        for bc in self.boundary_conditions:
            neumann_bc, dirichlet_bc = bc;
            self.deformer.clear();
            self.deformer.add_dirichlet_constraint(*dirichlet_bc);
            self.deformer.add_neumann_constraint(*neumann_bc[:2]);
            without_rigid_motion_constraint = len(dirichlet_bc[0]) != 0;
            displacement = self.deformer.solve(without_rigid_motion_constraint);
            self.displacements.append(displacement);

            stress = displacement_to_stress(self.assembler, displacement);
            stress = stress.reshape((-1, stress_size), order="C");
            stress_trace = np.sum(stress[:,0:dim], axis=1);
            self.stress_traces.append(stress_trace);

            applied_nodes = neumann_bc[0];
            applied_force = neumann_bc[1];
            applied_node_area = neumann_bc[2];

            force = np.zeros((self.mesh.num_vertices, dim));
            pressure = np.zeros((self.mesh.num_vertices, dim));
            if len(applied_nodes) > 0:
                force[applied_nodes] = applied_force;
                pressure[applied_nodes] = applied_force /\
                        np.array(applied_node_area)[:,np.newaxis];
            self.forces.append(force.ravel(order="C"));
            self.pressures.append(pressure.ravel(order="C"));

    @timethis
    def _initialize_coarse_mesh(self):
        self._create_coarse_mesh();
        self._create_coarse_mesh_assembler();
        self._extract_coarse_displacements();

    @timethis
    def _create_coarse_mesh(self):
        num_samples = 2;
        bbox_min, bbox_max = self.mesh.bbox;
        self.coarse_mesh = generate_box_mesh(bbox_min, bbox_max, num_samples);
        self.coarse_mesh.add_attribute("face_area");
        self.coarse_mesh.add_attribute("voxel_volume");

    @timethis
    def _create_coarse_mesh_assembler(self):
        coarse_material = Material(self.coarse_mesh.dim, None);
        coarse_assembler = PyAssembler.FEAssembler.create(
                self.coarse_mesh.raw_mesh, coarse_material.material);
        coarse_assembler = ElasticModel2.PyAssembler(
                self.coarse_mesh, coarse_assembler, coarse_material);
        self.coarse_assembler = coarse_assembler;

    @timethis
    def _extract_coarse_displacements(self):
        projector = AttributeProjection(self.coarse_mesh);
        self.coarse_displacements = projector.project(self.mesh, self.displacements);

    @timethis
    def _fit_material_parameters(self):
        raise NotImplementedError("This method is abstract");

    @property
    def mesh(self):
        return self.__mesh;

    @mesh.setter
    def mesh(self, mesh):
        self.__mesh = mesh;

    @property
    def bc_configs(self):
        raise NotImplementedError("This method is abstract");

