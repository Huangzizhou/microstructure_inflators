import numpy as np
from numpy.linalg import lstsq

import LinearElasticitySettings

from Mesh import Mesh
from mesh_io import load_mesh
import PyAssembler
from BoundaryCondition import BoundaryCondition
import ElasticModel2
from LinearElasticity import LinearElasticity
from ElasticityUtils import displacement_to_stress, total_energy,\
        force_to_pressure
from timethis import timethis

class OrthotropicFitter(object):
    def __init__(self, mesh):
        self.mesh = mesh;
        assembler = PyAssembler.FEAssembler.create_from_name(
                mesh.raw_mesh, "test_material");
        self.assembler = ElasticModel2.PyAssembler(mesh, assembler);
        self.deformer = LinearElasticity(self.mesh, self.assembler);
        self.displacements = [];
        self.stress_traces = [];
        self.forces = [];
        self.pressures = [];

    @timethis
    def fit(self):
        self.__generate_boundary_conditions();
        self.__deform_shape();
        self.__fit_orthotropic_parameters();

    @timethis
    def __generate_boundary_conditions(self):
        bd = BoundaryCondition(self.mesh);
        bc_configs = [
                self.__compression_x(),
                self.__compression_y(),
                self.__shear_x(),
                self.__compression_uniform() ];
        self.boundary_conditions = [bd.extract_from_dict(config)
                for config in bc_configs ];

    @timethis
    def __deform_shape(self):
        for bc in self.boundary_conditions:
            neumann_bc, dirichlet_bc = bc;
            self.deformer.clear();
            self.deformer.add_dirichlet_constraint(*dirichlet_bc);
            self.deformer.add_neumann_constraint(*neumann_bc);
            without_rigid_motion_constraint = len(dirichlet_bc[0]) != 0;
            displacement = self.deformer.solve(without_rigid_motion_constraint);
            self.displacements.append(displacement);

            applied_nodes = neumann_bc[0];
            applied_force = neumann_bc[1];
            force = np.zeros((self.mesh.num_vertices, self.mesh.dim));
            if len(applied_nodes) > 0:
                force[applied_nodes] = applied_force;
            self.forces.append(force.ravel(order="C"));
            pressure = force_to_pressure(self.mesh, force);
            self.pressures.append(pressure);

    @timethis
    def __fit_orthotropic_parameters(self):
        coeff = [];
        rhs = [];
        face_areas = self.mesh.get_attribute("face_area").ravel();
        for u in self.displacements:
            u = u.ravel();
            stress = displacement_to_stress(self.assembler, u);
            energy = total_energy(self.assembler, u);

            stress = stress.reshape((-1, 3), order="C");
            stress_trace = np.sum(stress[:,0:2], axis=1);
            energy = energy.ravel()[0];

            sigma_11 = stress[:,0].ravel();
            sigma_22 = stress[:,1].ravel();
            sigma_12 = stress[:,2].ravel();

            integrated_sigma_11_squared = np.multiply(sigma_11, sigma_11).dot(face_areas);
            integrated_sigma_22_squared = np.multiply(sigma_22, sigma_22).dot(face_areas);
            integrated_sigma_12_squared = np.multiply(sigma_12, sigma_12).dot(face_areas);
            integrated_sigma_11_time_22 = np.multiply(sigma_11, sigma_22).dot(face_areas);

            coeff.append([
                integrated_sigma_11_squared,
                integrated_sigma_22_squared,
                integrated_sigma_12_squared*2,
                integrated_sigma_11_time_22*2 ]);
            rhs.append(energy);
            self.stress_traces.append(stress_trace);

        parameter, residual, rank, singular_vals =\
                lstsq(coeff, rhs);
        self.orthotropic_parameter = parameter;
        self.residual_error = residual;
        self.condition_num = np.max(singular_vals) / np.min(singular_vals);
        if isinstance(self.residual_error, np.ndarray) and\
                len(self.residual_error) == 0:
            self.residual_error = 0.0;

    def __compression_x(self):
        eps = 1e-3;
        bbox_min, bbox_max = self.mesh.bbox;
        bc_config = {
                "neumann": [
                    {
                        "box": [bbox_min[0]-eps, bbox_min[0]+eps,
                            bbox_min[1], bbox_max[1]],
                        "force": [1.0, 0.0]
                    },
                    {
                        "box": [bbox_max[0]-eps, bbox_max[0]+eps,
                            bbox_min[1], bbox_max[1]],
                        "force": [-1.0, 0.0]
                    }
                    ]
                };
        return bc_config;

    def __compression_y(self):
        eps = 1e-3;
        bbox_min, bbox_max = self.mesh.bbox;
        bc_config = {
                "neumann": [
                    {
                        "box": [bbox_min[0], bbox_max[0],
                            bbox_min[1]-eps, bbox_min[1]+eps],
                        "force": [0.0, 1.0]
                    },
                    {
                        "box": [bbox_min[0], bbox_max[0],
                            bbox_max[1]-eps, bbox_max[1]+eps],
                        "force": [0.0, -1.0]
                    }
                    ]
                };
        return bc_config;

    def __shear_x(self):
        eps = 1e-3;
        bbox_min, bbox_max = self.mesh.bbox;
        bc_config = {
                "neumann": [
                    {
                        "box": [bbox_min[0], bbox_max[0],
                            bbox_max[1]-eps, bbox_max[1]+eps],
                        "force": [1.0, 0.0]
                    }
                ],
                "dirichlet": [
                    {
                        "box": [bbox_min[0], bbox_max[0],
                            bbox_min[1]-eps, bbox_min[1]+eps],
                        "fix": "all"
                    }
                    ]
                };
        return bc_config;

    def __compression_uniform(self):
        eps = 1e-3;
        bbox_min, bbox_max = self.mesh.bbox;
        bc_config = {
                "neumann": [
                    {
                        "box": [bbox_min[0]-eps, bbox_min[0]+eps,
                            bbox_min[1], bbox_max[1]],
                        "force": [1.0, 0.0]
                    },
                    {
                        "box": [bbox_max[0]-eps, bbox_max[0]+eps,
                            bbox_min[1], bbox_max[1]],
                        "force": [-1.0, 0.0]
                    },
                    {
                        "box": [bbox_min[0], bbox_max[0],
                            bbox_min[1]-eps, bbox_min[1]+eps],
                        "force": [0.0, 1.0]
                    },
                    {
                        "box": [bbox_min[0], bbox_max[0],
                            bbox_max[1]-eps, bbox_max[1]+eps],
                        "force": [0.0, -1.0]
                    }
                    ]
                };
        return bc_config;

    @property
    def mesh(self):
        return self.__mesh;

    @mesh.setter
    def mesh(self, mesh):
        if mesh.dim != 2:
            raise NotImplementedError("Only 2D mesh are supported for now.");
        self.__mesh = mesh;
        if not self.__mesh.has_attribute("face_area"):
            self.__mesh.add_attribute("face_area");

    @property
    def youngs_modulus(self):
        return [1.0 / self.orthotropic_parameter[0],
                1.0 / self.orthotropic_parameter[1] ];

    @property
    def poisson_ratio(self):
        young = self.youngs_modulus;
        return [self.orthotropic_parameter[3]*young[0],
                self.orthotropic_parameter[3]*young[1] ];

    @property
    def shear_modulus(self):
        return 1.0 / self.orthotropic_parameter[2];

