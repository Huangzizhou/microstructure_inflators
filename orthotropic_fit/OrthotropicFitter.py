import numpy as np
from numpy.linalg import lstsq, inv, norm

import LinearElasticitySettings
import PyMeshUtils

from Mesh import Mesh
from mesh_io import load_mesh, form_mesh
import PyAssembler
from BoundaryConditionExtractor import BoundaryConditionExtractor
import ElasticModel2
from BoxMeshGenerator import generate_box_mesh
from ElasticityUtils import displacement_to_stress, displacement_to_strain,\
        total_energy
from LinearElasticity import LinearElasticity
from Material import Material
import PredefinedBoundaryConditions as PredefinedBC
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
        bbox_min, bbox_max = self.mesh.bbox;
        if self.mesh.dim == 2:
            eps = 1e-3;
            self.bc_configs = [
                    PredefinedBC.compress_x(bbox_min, bbox_max, eps),
                    PredefinedBC.compress_y(bbox_min, bbox_max, eps),
                    PredefinedBC.compress_xy(bbox_min, bbox_max, eps),
                    ];
        elif self.mesh.dim == 3:
            eps = norm(bbox_max - bbox_min) * 0.01;
            self.bc_configs = [
                    PredefinedBC.compress_x(bbox_min, bbox_max, eps),
                    PredefinedBC.compress_y(bbox_min, bbox_max, eps),
                    PredefinedBC.compress_z(bbox_min, bbox_max, eps),
                    PredefinedBC.compress_xy(bbox_min, bbox_max, eps),
                    PredefinedBC.compress_yz(bbox_min, bbox_max, eps),
                    PredefinedBC.compress_zx(bbox_min, bbox_max, eps),
                    ];
        else:
            raise RuntimeError("Unsupported dim: {}".format(self.mesh.dim));

        self.boundary_conditions = [];
        for config in self.bc_configs:
            bd.clear();
            bd.extract_from_dict(config);
            self.boundary_conditions.append([
                bd.neumann_bc, bd.dirichlet_bc]);

    @timethis
    def __deform_shape(self):
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

    def __fit_orthotropic_parameters(self):
        if self.mesh.dim == 2:
            return self.__fit_orthotropic_parameters_2D();
        elif self.mesh.dim == 3:
            return self.__fit_orthotropic_parameters_3D();
        else:
            raise RuntimeError("Unsupported dim: {}".format(self.mesh.dim));

    @timethis
    def __fit_orthotropic_parameters_2D(self):
        dim = self.mesh.dim;
        tensor_size = dim*(dim+1)/2;
        assert(tensor_size == 3);
        bbox_min, bbox_max = self.mesh.bbox;
        num_samples = 2;
        eps = 1e-6;
        self.coarse_mesh = generate_box_mesh(bbox_min-eps, bbox_max+eps, num_samples);
        self.coarse_mesh.add_attribute("face_area");
        face_areas = self.coarse_mesh.get_attribute("face_area").ravel();

        self.coarse_displacements = self.__extract_coarse_displacements();
        coarse_material = Material(self.coarse_mesh.dim, None);
        coarse_assembler = PyAssembler.FEAssembler.create(
                self.coarse_mesh.raw_mesh, coarse_material.material);
        coarse_assembler = ElasticModel2.PyAssembler(
                self.coarse_mesh, coarse_assembler, coarse_material);
        K = self.assembler.stiffness_matrix;

        coeff = [];
        rhs = [];
        for u1,coarse_u1 in zip(self.displacements, self.coarse_displacements):
            strain_1 = displacement_to_strain(coarse_assembler, coarse_u1);
            strain_1 = strain_1.reshape((-1, tensor_size), order="C");

            strain_1_xx = strain_1[:,0].ravel();
            strain_1_yy = strain_1[:,1].ravel();
            strain_1_xy = strain_1[:,2].ravel();

            for u2,coarse_u2 in zip(self.displacements, self.coarse_displacements):
                energy = np.dot(u1, K*u2);

                strain_2 = displacement_to_strain(coarse_assembler, coarse_u2);
                strain_2 = strain_2.reshape((-1, tensor_size), order="C");

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
        self.orthotropic_parameter = np.array(parameter);
        self.residual_error = residual;
        self.condition_num = np.max(singular_vals) / np.min(singular_vals);
        if isinstance(self.residual_error, np.ndarray):
            if len(self.residual_error) == 0:
                self.residual_error = 0.0;
            else:
                self.residual_error = np.max(self.residual_error);

    @timethis
    def __fit_orthotropic_parameters_3D(self):
        dim = self.mesh.dim;
        tensor_size = dim*(dim+1)/2;
        bbox_min, bbox_max = self.mesh.bbox;
        num_samples = 3;
        eps = 1e-6;
        self.coarse_mesh = generate_box_mesh(bbox_min-eps, bbox_max+eps, num_samples);
        self.coarse_mesh.add_attribute("voxel_volume");
        voxel_volumes = self.coarse_mesh.get_attribute("voxel_volume").ravel();

        self.coarse_displacements = self.__extract_coarse_displacements();
        coarse_material = Material(self.coarse_mesh.dim, None);
        coarse_assembler = PyAssembler.FEAssembler.create(
                self.coarse_mesh.raw_mesh, coarse_material.material);
        coarse_assembler = ElasticModel2.PyAssembler(
                self.coarse_mesh, coarse_assembler, coarse_material);
        K = self.assembler.stiffness_matrix;

        coeff = [];
        rhs = [];
        for u1,coarse_u1 in zip(self.displacements, self.coarse_displacements):
            strain_1 = displacement_to_strain(coarse_assembler, coarse_u1);
            strain_1 = strain_1.reshape((-1, tensor_size), order="C");

            strain_1_xx = strain_1[:,0].ravel();
            strain_1_yy = strain_1[:,1].ravel();
            strain_1_zz = strain_1[:,2].ravel();
            strain_1_xy = strain_1[:,3].ravel();
            strain_1_xz = strain_1[:,4].ravel();
            strain_1_yz = strain_1[:,5].ravel();
            for u2,coarse_u2 in zip(self.displacements, self.coarse_displacements):
                energy = np.dot(u1, K*u2);

                strain_2 = displacement_to_strain(coarse_assembler, coarse_u2);
                strain_2 = strain_2.reshape((-1, tensor_size), order="C");

                strain_2_xx = strain_2[:,0].ravel();
                strain_2_yy = strain_2[:,1].ravel();
                strain_2_zz = strain_2[:,2].ravel();
                strain_2_xy = strain_2[:,3].ravel();
                strain_2_xz = strain_2[:,4].ravel();
                strain_2_yz = strain_2[:,5].ravel();

                coeff_11 = np.multiply(strain_1_xx, strain_2_xx).dot(voxel_volumes);
                coeff_22 = np.multiply(strain_1_yy, strain_2_yy).dot(voxel_volumes);
                coeff_33 = np.multiply(strain_1_zz, strain_2_zz).dot(voxel_volumes);
                coeff_44 = np.multiply(strain_1_xy, strain_2_xy).dot(voxel_volumes) * 2;
                coeff_55 = np.multiply(strain_1_xz, strain_2_xz).dot(voxel_volumes) * 2;
                coeff_66 = np.multiply(strain_1_yz, strain_2_yz).dot(voxel_volumes) * 2;
                coeff_12 = (np.multiply(strain_1_xx, strain_2_yy) +\
                        np.multiply(strain_1_yy, strain_2_xx)).dot(voxel_volumes)
                coeff_13 = (np.multiply(strain_1_xx, strain_2_zz) +\
                        np.multiply(strain_1_zz, strain_2_xx)).dot(voxel_volumes)
                coeff_23 = (np.multiply(strain_1_yy, strain_2_zz) +\
                        np.multiply(strain_1_zz, strain_2_yy)).dot(voxel_volumes)

                coeff.append([
                    coeff_11, coeff_22, coeff_33,
                    coeff_44, coeff_55, coeff_66,
                    coeff_12, coeff_13, coeff_23]);
                rhs.append(energy);

        parameter, residual, rank, singular_vals =\
                lstsq(coeff, rhs);
        C = np.array([
            [parameter[0], parameter[6], parameter[7],          0.0,          0.0,          0.0],
            [parameter[6], parameter[1], parameter[8],          0.0,          0.0,          0.0],
            [parameter[7], parameter[8], parameter[2],          0.0,          0.0,          0.0],
            [         0.0,          0.0,          0.0, parameter[3],          0.0,          0.0],
            [         0.0,          0.0,          0.0,          0.0, parameter[4],          0.0],
            [         0.0,          0.0,          0.0,          0.0,          0.0, parameter[5]]]);

        S = inv(C);
        parameter = [
                S[0,0], S[1,1], S[2,2],
                S[3,3], S[4,4], S[5,5],
                S[0,1], S[0,2], S[1,2]];
        self.orthotropic_parameter = np.array(parameter);
        self.residual_error = residual;
        self.condition_num = np.max(singular_vals) / np.min(singular_vals);
        if isinstance(self.residual_error, np.ndarray):
            if len(self.residual_error) == 0:
                self.residual_error = 0.0;
            else:
                self.residual_error = np.max(self.residual_error);

    @timethis
    def __extract_coarse_displacements_old(self, mesh):
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
    def __extract_coarse_displacements(self):
        dim = self.mesh.dim;
        if dim == 2:
            coarse_elems = self.coarse_mesh.faces;
            vertex_per_elem = self.coarse_mesh.vertex_per_face;
        elif dim == 3:
            coarse_elems = self.coarse_mesh.voxels;
            vertex_per_elem = self.coarse_mesh.vertex_per_voxel;
        else:
            raise NotImplementedError("Only 2D and 3D are supported.");

        point_locator = PyMeshUtils.PointLocator(self.coarse_mesh.raw_mesh);
        point_locator.locate(self.mesh.vertices);
        element_indices = point_locator.get_enclosing_voxels().ravel();
        barycentric_coords = point_locator.get_barycentric_coords();

        barycentric_matrix = np.zeros((
            self.mesh.num_vertices, self.coarse_mesh.num_vertices));
        for i in range(self.mesh.num_vertices):
            barycentric_matrix[i, coarse_elems[element_indices[i]]] =\
                    barycentric_coords[i];

        coarse_displacements = [];
        for u in self.displacements:
            u = u.reshape((-1, dim), order="C");
            coarse_u, residual, rank, singular_vals  = lstsq(barycentric_matrix, u);
            coarse_displacements.append(coarse_u.ravel(order="C"));
        return coarse_displacements;

    @timethis
    def __extract_coarse_displacements_old(self):
        dim = self.mesh.dim;
        coarse_vertices = self.coarse_mesh.vertices;
        if dim == 2:
            coarse_elems = self.coarse_mesh.faces;
            vertex_per_elem = self.coarse_mesh.vertex_per_face;
        elif dim == 3:
            coarse_elems = self.coarse_mesh.voxels;
            vertex_per_elem = self.coarse_mesh.vertex_per_voxel;
        else:
            raise NotImplementedError("Only 2D and 3D are supported.");

        assert(coarse_elems.shape[1] == vertex_per_elem);

        barycentric_solver = [];
        for elem in coarse_elems:
            T = np.array(coarse_vertices[elem[:dim]] -
                    coarse_vertices[elem[-1]]).T;
            barycentric_solver.append(inv(T));
        last_vertex = coarse_vertices[coarse_elems[:,-1]];
        assert(len(barycentric_solver) == len(coarse_elems));

        fine_vertices = self.mesh.vertices;
        barycentric_matrix = np.zeros((
            self.mesh.num_vertices, self.coarse_mesh.num_vertices));
        for i,v in enumerate(fine_vertices):
            barycentric_coords = np.array([np.sum(M * (v-lv), axis=1)
                    for M,lv in zip(barycentric_solver, last_vertex)]);
            last_coord = 1.0 - np.sum(barycentric_coords, axis=1).reshape((-1,1));
            barycentric_coords = np.hstack((barycentric_coords, last_coord));

            for coord, elem in zip (barycentric_coords, coarse_elems):
                interpolated_v = np.sum(coarse_vertices[elem].T * coord,
                        axis=1).ravel();
                assert(norm(v-interpolated_v) < 1e-3);

            inside = np.all(np.logical_and(
                    barycentric_coords >= 0.0,
                    barycentric_coords <= 1.0), axis=1);
            assert(np.any(inside));
            row = barycentric_matrix[i];
            row[coarse_elems[inside]] = barycentric_coords[inside];
            barycentric_matrix[i] = row.ravel(order="C");

        coarse_displacements = [];
        for u in self.displacements:
            u = u.reshape((-1, dim), order="C");
            coarse_u, residual, rank, singular_vals  = lstsq(barycentric_matrix, u);
            coarse_displacements.append(coarse_u.ravel(order="C"));
        return coarse_displacements;



    @property
    def mesh(self):
        return self.__mesh;

    @mesh.setter
    def mesh(self, mesh):
        self.__mesh = mesh;

    @property
    def youngs_modulus(self):
        dim = self.mesh.dim;
        return 1.0 / self.orthotropic_parameter[:dim];

    @property
    def poisson_ratio(self):
        dim = self.mesh.dim;
        young = self.youngs_modulus;
        if dim == 2:
            return np.array([
                -self.orthotropic_parameter[3] * young[0], # v_xy
                -self.orthotropic_parameter[3] * young[1]  # v_yx
                ]);
        elif dim == 3:
            return np.array( [
                -self.orthotropic_parameter[6] * young[0], # v_xy
                -self.orthotropic_parameter[6] * young[1], # v_yx
                -self.orthotropic_parameter[7] * young[0], # v_xz
                -self.orthotropic_parameter[7] * young[2], # v_zx
                -self.orthotropic_parameter[8] * young[1], # v_yz
                -self.orthotropic_parameter[8] * young[2], # v_zy
                ]);

    @property
    def shear_modulus(self):
        dim = self.mesh.dim;
        if dim == 2:
            shear = 0.5 / self.orthotropic_parameter[2:3];
        elif dim == 3:
            shear = 0.5 / self.orthotropic_parameter[3:6];
        else:
            raise RuntimeError("Unsupported dim: {}".format(self.mesh.dim));

        if np.any(shear <= 0):
            import warnings
            warnings.warn("Negative shear modulus: {}".format(shear));
        return abs(shear);

