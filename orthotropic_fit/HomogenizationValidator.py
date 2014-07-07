import LinearElasticitySettings

import numpy as np

from mesh_io import form_mesh

from BoundaryConditionExtractor import BoundaryConditionExtractor
from ElasticityUtils import displacement_to_stress
import ElasticModel2
from BoxMeshGenerator import generate_box_mesh
from LinearElasticity import LinearElasticity
import PyAssembler
from timethis import timethis

class HomogenizationValidator(object):
    def __init__(self, mesh, material):
        self.input_mesh = mesh;
        self.material = material;
        self.__generate_bbox_mesh();
        self.__initialize_assembler();
        self.__initialize_linear_elasticity();
        self.displacements = [];
        self.stress_traces = [];
        self.forces = [];
        self.pressures = [];

    @timethis
    def simulate(self, bc_configs):
        dim = self.mesh.dim;
        tensor_size = dim * (dim+1) / 2;
        bd = BoundaryConditionExtractor(self.mesh);
        boundary_conditions = [];
        for config in bc_configs:
            bd.clear();
            bd.extract_from_dict(config);
            neumann_bc = bd.neumann_bc;
            dirichlet_bc = bd.dirichlet_bc;

            self.deformer.clear();
            self.deformer.add_dirichlet_constraint(*dirichlet_bc[:2]);
            self.deformer.add_neumann_constraint(*neumann_bc[:2]);
            without_rigid_motion_constraint = len(dirichlet_bc[0]) != 0;
            displacement = self.deformer.solve(without_rigid_motion_constraint);
            self.displacements.append(displacement);

            stress = displacement_to_stress(self.assembler, displacement);
            stress = stress.reshape((-1, tensor_size), order="C");
            stress_trace = np.sum(stress[:,0:dim], axis=1);
            self.stress_traces.append(stress_trace);

            applied_nodes = neumann_bc[0];
            applied_force = neumann_bc[1];
            applied_node_area = neumann_bc[2];
            force = np.zeros((self.mesh.num_vertices, self.mesh.dim));
            pressure = np.zeros((self.mesh.num_vertices, self.mesh.dim));
            if len(applied_nodes) > 0:
                force[applied_nodes] = applied_force;
                pressure[applied_nodes] = applied_force /\
                        np.array(applied_node_area)[:,np.newaxis];
            self.forces.append(force.ravel(order="C"));
            self.pressures.append(pressure.ravel(order="C"));

    @timethis
    def __generate_bbox_mesh(self):
        bbox_min, bbox_max = self.input_mesh.bbox;
        num_samples = 10;
        self.mesh = generate_box_mesh(bbox_min, bbox_max, num_samples);

    @timethis
    def __initialize_assembler(self):
        assembler = PyAssembler.FEAssembler.create(self.mesh.raw_mesh, self.material);
        self.assembler = ElasticModel2.PyAssembler(self.mesh, assembler,
                self.material);

    @timethis
    def __initialize_linear_elasticity(self):
        self.deformer = LinearElasticity(self.mesh, self.assembler);

