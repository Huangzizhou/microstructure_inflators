import numpy as np
from numpy.linalg import lstsq, inv, norm

import LinearElasticitySettings

from Mesh import Mesh
from mesh_io import load_mesh, form_mesh
import PyAssembler
from BoundaryConditionExtractor import BoundaryConditionExtractor
import ElasticModel2
from generate_box_mesh import generate_box_mesh
from ElasticityUtils import displacement_to_stress, displacement_to_strain,\
        total_energy, force_to_pressure
from LinearElasticity import LinearElasticity
from Material import Material
from timethis import timethis

class OrthotropicFitter(object):
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
        self.__generate_boundary_conditions();
        self.__deform_shape();
        self.__fit_orthotropic_parameters();

    @timethis
    def __generate_boundary_conditions(self):
        bd = BoundaryConditionExtractor(self.mesh);
        self.bc_configs = [
                self.__compression_x(),
                self.__compression_y(),
                self.__compression_xy(),
                #self.__shear_x(),
                #self.__shear_y(),
                #self.__compression_uniform()
                ];
        self.boundary_conditions = [];
        for config in self.bc_configs:
            bd.clear();
            bd.extract_from_dict(config);
            self.boundary_conditions.append([
                bd.neumann_bc, bd.dirichlet_bc]);

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
    def __fit_orthotropic_parameters(self):
        bbox_min, bbox_max = self.mesh.bbox;
        num_samples = 2;
        mesh = generate_box_mesh(bbox_min, bbox_max, num_samples);
        mesh.add_attribute("face_area");
        face_areas = mesh.get_attribute("face_area").ravel();

        coarse_displacements = self.__extract_coarse_displacements(mesh);
        coarse_material = Material(mesh.dim, None);
        coarse_assembler = PyAssembler.FEAssembler.create(
                mesh.raw_mesh, coarse_material.material);
        coarse_assembler = ElasticModel2.PyAssembler(
                mesh, coarse_assembler, coarse_material);
        K = self.assembler.stiffness_matrix;

        coeff = [];
        rhs = [];
        for u1,coarse_u1 in zip(self.displacements, coarse_displacements):
            for u2,coarse_u2 in zip(self.displacements, coarse_displacements):
                energy = np.dot(u1, K*u2);
                strain_1 = displacement_to_strain(coarse_assembler, coarse_u1);
                strain_2 = displacement_to_strain(coarse_assembler, coarse_u2);
                strain_1 = strain_1.reshape((-1, 3), order="C");
                strain_2 = strain_2.reshape((-1, 3), order="C");

                strain_1_xx = strain_1[:,0].ravel();
                strain_1_yy = strain_1[:,1].ravel();
                strain_1_xy = strain_1[:,2].ravel();

                strain_2_xx = strain_2[:,0].ravel();
                strain_2_yy = strain_2[:,1].ravel();
                strain_2_xy = strain_2[:,2].ravel();

                coeff_11 = np.multiply(strain_1_xx, strain_2_xx).dot(face_areas);
                coeff_22 = np.multiply(strain_1_yy, strain_2_yy).dot(face_areas);
                coeff_33 = np.multiply(strain_1_xy, strain_2_xy).dot(face_areas) * 2;
                coeff_12 = (np.multiply(strain_1_xx, strain_2_yy) + \
                        np.multiply(strain_1_yy, strain_2_xx)).dot(face_areas);

                coeff.append([coeff_11, coeff_22, coeff_33, coeff_12]);
                rhs.append(energy);

        parameter, residual, rank, singular_vals =\
                lstsq(coeff, rhs);
        C = np.array([
            [parameter[0], parameter[3],          0.0],
            [parameter[3], parameter[1],          0.0],
            [         0.0,          0.0, parameter[2]] ]);
        S = inv(C);
        parameter = [S[0,0], S[1,1], S[2,2], S[0,1]];
        self.orthotropic_parameter = parameter;
        self.residual_error = residual;
        self.condition_num = np.max(singular_vals) / np.min(singular_vals);
        if isinstance(self.residual_error, np.ndarray):
            if len(self.residual_error) == 0:
                self.residual_error = 0.0;
            else:
                self.residual_error = np.max(self.residual_error);

    @timethis
    def __extract_coarse_displacements(self, mesh):
        vertices = self.mesh.vertices;

        nearest_idx = np.zeros(mesh.num_vertices, dtype=int);
        for i,v in enumerate(mesh.vertices):
            distances = norm(vertices - v, axis=1);
            nearest_idx[i] = np.argmin(distances);

        coarse_displacements = [];
        for u in self.displacements:
            u = u.reshape((-1, mesh.dim), order="C");
            coarse_u = u[nearest_idx, :].ravel(order="C");
            coarse_displacements.append(coarse_u);

        return coarse_displacements;

    @timethis
    def __fit_orthotropic_parameters_old(self):
        face_areas = self.mesh.get_attribute("face_area").ravel();
        total_area = np.sum(face_areas);
        bbox = self.mesh.bbox;
        bbox_area = np.prod(bbox[1] - bbox[0]);
        K = self.assembler.stiffness_matrix;

        coeff = [];
        rhs = [];
        for u1 in self.displacements:
            for u2 in self.displacements:
                energy = np.dot(u1, K*u2);
                strain_1 = displacement_to_strain(self.assembler, u1);
                strain_2 = displacement_to_strain(self.assembler, u2);
                strain_1 = strain_1.reshape((-1, 3), order="C");
                strain_2 = strain_2.reshape((-1, 3), order="C");

                strain_1_xx = strain_1[:,0].ravel();
                strain_1_yy = strain_1[:,1].ravel();
                strain_1_xy = strain_1[:,2].ravel();

                strain_2_xx = strain_2[:,0].ravel();
                strain_2_yy = strain_2[:,1].ravel();
                strain_2_xy = strain_2[:,2].ravel();

                coeff_11 = np.multiply(strain_1_xx, strain_2_xx).dot(face_areas);
                coeff_22 = np.multiply(strain_1_yy, strain_2_yy).dot(face_areas);
                coeff_33 = np.multiply(strain_1_xy, strain_2_xy).dot(face_areas) * 2;
                coeff_12 = (np.multiply(strain_1_xx, strain_2_yy) + \
                        np.multiply(strain_1_yy, strain_2_xx)).dot(face_areas);

                coeff_11 *= bbox_area / total_area;
                coeff_22 *= bbox_area / total_area;
                coeff_33 *= bbox_area / total_area;
                coeff_12 *= bbox_area / total_area;

                coeff.append([coeff_11, coeff_22, coeff_33, coeff_12]);
                rhs.append(energy);

        parameter, residual, rank, singular_vals =\
                lstsq(coeff, rhs);
        C = np.array([
            [parameter[0], parameter[3],          0.0],
            [parameter[3], parameter[1],          0.0],
            [         0.0,          0.0, parameter[2]] ]);
        S = inv(C);
        parameter = [S[0,0], S[1,1], S[2,2], S[0,1]];
        self.orthotropic_parameter = parameter;
        self.residual_error = residual;
        self.condition_num = np.max(singular_vals) / np.min(singular_vals);
        if isinstance(self.residual_error, np.ndarray):
            if len(self.residual_error) == 0:
                self.residual_error = 0.0;
            else:
                self.residual_error = np.max(self.residual_error);

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

    def __shear_y(self):
        eps = 1e-3;
        bbox_min, bbox_max = self.mesh.bbox;
        bc_config = {
                "neumann": [
                    {
                        "box": [bbox_max[0]-eps, bbox_max[0]+eps,
                            bbox_min[1], bbox_max[1]],
                        "force": [0.0, 1.0]
                    }
                ],
                "dirichlet": [
                    {
                        "box": [bbox_min[0]-eps, bbox_min[0]+eps,
                            bbox_min[1], bbox_max[1]],
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

    def __compression_xy(self):
        eps = 1e-3;
        bbox_min, bbox_max = self.mesh.bbox;
        bc_config = {
                "neumann": [
                    {
                        "box": [bbox_min[0]-eps, bbox_min[0]+eps,
                            bbox_min[1], bbox_max[1]],
                        "force": [0.0, -1.0]
                    },
                    {
                        "box": [bbox_max[0]-eps, bbox_max[0]+eps,
                            bbox_min[1], bbox_max[1]],
                        "force": [0.0, 1.0]
                    },
                    {
                        "box": [bbox_min[0], bbox_max[0],
                            bbox_min[1]-eps, bbox_min[1]+eps],
                        "force": [-1.0, 0.0]
                    },
                    {
                        "box": [bbox_min[0], bbox_max[0],
                            bbox_max[1]-eps, bbox_max[1]+eps],
                        "force": [1.0, 0.0]
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
        return np.array([1.0 / self.orthotropic_parameter[0],
            1.0 / self.orthotropic_parameter[1] ]);

    @property
    def poisson_ratio(self):
        young = self.youngs_modulus;
        return np.array([-self.orthotropic_parameter[3]*young[0],
                -self.orthotropic_parameter[3]*young[1] ]);

    @property
    def shear_modulus(self):
        shear = np.array([0.5 / self.orthotropic_parameter[2]]);
        if shear <= 0:
            import warnings
            warnings.warn("Negative shear modulus: {}".format(shear));
        return abs(shear);

