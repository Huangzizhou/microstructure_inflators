import LinearElasticitySettings

import numpy as np

from mesh_io import form_mesh

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
        for config in bc_configs:
            self.deformer.clear();
            self.deformer.add_boundary_condition_from_dict(config);
            without_rigid_motion_constraint = len(self.deformer.fixed_nodes) > 0;
            displacement = self.deformer.solve(without_rigid_motion_constraint);
            self.displacements.append(displacement);

            stress = displacement_to_stress(self.assembler, displacement);
            stress = stress.reshape((-1, tensor_size), order="C");
            stress_trace = np.sum(stress[:,0:dim], axis=1);
            self.stress_traces.append(stress_trace);

            force = self.deformer.force;
            pressure = self.deformer.traction;
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

