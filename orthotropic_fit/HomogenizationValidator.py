import LinearElasticitySettings

import numpy as np

from mesh_io import form_mesh

from BoundaryCondition import BoundaryCondition
from ElasticityUtils import displacement_to_stress, force_to_pressure
import ElasticModel2
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
        bd = BoundaryCondition(self.mesh);
        boundary_conditions = [bd.extract_from_dict(config)
                for config in bc_configs];
        for bc in boundary_conditions:
            neumann_bc, dirichlet_bc = bc;
            self.deformer.clear();
            self.deformer.add_dirichlet_constraint(*dirichlet_bc);
            self.deformer.add_neumann_constraint(*neumann_bc);
            without_rigid_motion_constraint = len(dirichlet_bc[0]) != 0;
            displacement = self.deformer.solve(without_rigid_motion_constraint);
            self.displacements.append(displacement);

            stress = displacement_to_stress(self.assembler, displacement);
            stress = stress.reshape((-1, 3), order="C");
            stress_trace = np.sum(stress[:,0:2], axis=1);
            self.stress_traces.append(stress_trace);

            applied_nodes = neumann_bc[0];
            applied_force = neumann_bc[1];
            force = np.zeros((self.mesh.num_vertices, self.mesh.dim));
            if len(applied_nodes) > 0:
                force[applied_nodes] = applied_force;
            self.forces.append(force.ravel(order="C"));
            pressure = force_to_pressure(self.mesh, force);
            self.pressures.append(pressure);

    @timethis
    def __generate_bbox_mesh(self):
        if self.input_mesh.dim != 2:
            raise NotImplementedError("Only 2D mesh are supported for now.");
        bbox_min, bbox_max = self.input_mesh.bbox;
        num_samples = 100;
        vertices = [];
        faces = [];
        for i in range(num_samples+1):
            for j in range(num_samples+1):
                x_ratio = float(i) / float(num_samples);
                y_ratio = float(j) / float(num_samples);
                x = x_ratio * bbox_max[0] + (1.0 - x_ratio) * bbox_min[0];
                y = y_ratio * bbox_max[1] + (1.0 - y_ratio) * bbox_min[1];
                vertices.append([x,y]);

        row_size = num_samples+1;
        for i in range(num_samples):
            for j in range(num_samples):
                idx = [i*row_size+j,
                        i*row_size+j+1,
                        (i+1)*row_size+j,
                        (i+1)*row_size+j+1 ];
                faces.append([idx[0], idx[2], idx[3]]);
                faces.append([idx[0], idx[3], idx[1]]);

        vertices = np.array(vertices, dtype=float);
        faces = np.array(faces, dtype=int);
        self.mesh = form_mesh(vertices, faces);

    @timethis
    def __initialize_assembler(self):
        assembler = PyAssembler.FEAssembler.create(self.mesh.raw_mesh, self.material);
        self.assembler = ElasticModel2.PyAssembler(self.mesh, assembler);

    @timethis
    def __initialize_linear_elasticity(self):
        self.deformer = LinearElasticity(self.mesh, self.assembler);

