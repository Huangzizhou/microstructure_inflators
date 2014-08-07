import hashlib
import numpy as np
from numpy.linalg import norm
import scipy.sparse
import scipy.sparse.linalg

import LinearElasticitySettings
import PyAssembler
from timethis import timethis
from MatrixUtils import format

from OptimizationSetting import OptimizationSetting

class IsotropicMatOptSetting(OptimizationSetting):
    def __init__(self, mesh, bc_extractor, init_young, init_poisson):
        self.mesh = mesh;

        self.__initialize_optimization_parameters();
        self.__initialize_neumann_bc(bc_extractor);
        self.__initialize_dirichlet_bc(bc_extractor);
        self.__initialize_material_parameters(init_young, init_poisson);
        self.__initialize_elasticity_gradient();
        self.__initialize_materials();
        self.__initialize_assembler();
        self.__initialize_matrices();
        self.__initialize_rigid_motion_rhs();
        self.__initialize_history();

    def evaluate(self, parameters):
        obj, grad = self.__evaluate_objective_and_gradient(parameters);
        return obj, grad;

    def log_iteration(self, parameters):
        key = self.__parameters_to_key(parameters);
        if key in self.cache:
            idx = self.cache[key];
            self.iteration_indices.append(idx);
            print("iteration: {}  obj={}  grad={}".format(
                len(self.iteration_indices),
                self.objective_history[idx],
                norm(self.gradient_history[idx])));
        else:
            assert(False);

    def __initialize_optimization_parameters(self):
        self.num_dof = self.mesh.num_vertices * self.mesh.dim;
        self.source_term = np.zeros(self.num_dof);
        self.target_displacement = np.zeros(self.num_dof);

    def __initialize_neumann_bc(self, bc_extractor):
        self.applied_idx = bc_extractor.applied_node_idx;
        self.applied_traction = np.multiply(
                bc_extractor.applied_node_traction, 
                np.array(bc_extractor.applied_node_area)[:,np.newaxis]);

        assert(self.mesh.dim == self.applied_traction.shape[1]);

        neumann_term = np.zeros((self.mesh.num_vertices, self.mesh.dim));
        neumann_term[self.applied_idx] = self.applied_traction;
        self.source_term += neumann_term.ravel(order="C");

    def __initialize_dirichlet_bc(self, bc_extractor):
        self.fixed_idx, self.fixed_position, self.fixed_node_area\
                = bc_extractor.dirichlet_bc;

        dirichlet_term = np.zeros((self.mesh.num_vertices, self.mesh.dim));
        dirichlet_term[self.fixed_idx] = self.fixed_position;
        self.target_displacement += dirichlet_term.ravel(order="C");
        self.target_areas = np.zeros(self.mesh.num_vertices);
        self.target_areas[self.fixed_idx] = self.fixed_node_area;

    def __initialize_elasticity_gradient(self):
        if self.mesh.dim == 2:
            self.C_lambda = np.array([
                [1, 1, 0],
                [1, 1, 0],
                [0, 0, 0] ], dtype=float);
            self.C_mu = np.diag([2.0, 2.0, 1.0]);

            self.C_E = lambda E, v: \
                    self.C_lambda * v / (1 - v**2) +\
                    self.C_mu * 1.0 / (2 + 2 * v);
            self.C_v = lambda E, v: \
                    self.C_lambda * (E / (1-v**2) + 2*E*v**2/(1-v**2)**2) -\
                    self.C_mu * 2*E / (2+2*v)**2;

            self.grad_young_field_name = "C_grad_young";
            self.mesh.add_attribute(self.grad_young_field_name);
            self.grad_poisson_field_name = "C_grad_poisson";
            self.mesh.add_attribute(self.grad_poisson_field_name);

            self.__update_elasticity_gradient();
        elif self.mesh.dim == 3:
            # TODO
            raise NotImplementedError("3D is not yet supported");
        else:
            raise NotImplementedError("{}D is not yet supported"\
                    .format(self.mesh.dim));

    def __initialize_material_parameters(self, init_young, init_poisson):
        self.density = 1.0

        self.young_field_name = "young";
        self.mesh.add_attribute(self.young_field_name);
        self.mesh.set_attribute(self.young_field_name, init_young);

        self.poisson_field_name = "poisson";
        self.mesh.add_attribute(self.poisson_field_name);
        self.mesh.set_attribute(self.poisson_field_name, init_poisson);

    def __initialize_materials(self):
        self.hetero_material =\
                PyAssembler.Material.create_element_wise_isotropic(
                        self.density,
                        self.mesh.raw_mesh, 
                        self.young_field_name,
                        self.poisson_field_name);

        self.grad_young_material = \
                PyAssembler.Material.create_element_wise_symmetric(
                        self.density, self.mesh.raw_mesh,
                        self.grad_young_field_name);

        self.grad_poisson_material = \
                PyAssembler.Material.create_element_wise_symmetric(
                        self.density, self.mesh.raw_mesh,
                        self.grad_poisson_field_name);

    def __initialize_assembler(self):
        self.assembler = PyAssembler.FEAssembler.create(
                self.mesh.raw_mesh, self.hetero_material);

    def __initialize_matrices(self):
        dim = self.mesh.dim;
        self.stiffness = format(self.assembler.assemble("stiffness"));
        self.rigid_motion = format(self.assembler.assemble("rigid_motion"));
        self.displacement_strain = format(self.assembler.assemble("displacement_strain"));

        flattened_tensor_size = dim*(dim+1)/2;
        size = self.mesh.num_elements * flattened_tensor_size;

        diagonal = np.ones((self.mesh.num_elements, flattened_tensor_size));
        diagonal[:,dim:] = 2;
        diagonal = diagonal.ravel(order="C");
        self.strain_doubler = scipy.sparse.dia_matrix(
                (diagonal, 0), shape=(size, size));

        volumes = self.mesh.element_volumes;
        self.element_volumes = scipy.sparse.dia_matrix(
                (np.repeat(volumes, flattened_tensor_size), 0),
                shape=(size, size));

        num_faces = self.mesh.num_faces;
        i_idx = [];
        j_idx = [];
        val = [];
        self.mesh.raw_mesh.enable_face_connectivity();
        for i in range(num_faces):
            adj_faces = self.mesh.raw_mesh.get_face_adjacent_faces(i).ravel();
            for j in adj_faces:
                i_idx.append(i);
                j_idx.append(i);
                val.append(1.0);

                i_idx.append(i);
                j_idx.append(j);
                val.append(-1.0);

        self.element_connectivity = scipy.sparse.coo_matrix((val, (i_idx,
            j_idx)), shape=(num_faces, num_faces)).tocsc() * 1e-5;

    def __initialize_rigid_motion_rhs(self):
        dim = self.mesh.dim;
        K = self.stiffness;
        i = []; j = []; val = []; rhs = [];
        for counter, idx, pos in zip(
                range(len(self.fixed_idx)),
                self.fixed_idx,
                self.fixed_position.reshape((-1, dim))):
            i += range(counter*dim, (counter+1)*dim);
            j += range(idx*dim, (idx+1)*dim)
            val += [1.0] * dim;
            rhs += pos.tolist();
        R = scipy.sparse.coo_matrix((val, (i, j)),
                shape=(len(self.fixed_idx)*dim, self.num_dof)).tocsc();
        A = scipy.sparse.bmat([[K, R.T], [R, None]]).tocsc();
        b = np.zeros(A.shape[0]);
        b[self.num_dof:] = rhs;

        sol = scipy.sparse.linalg.spsolve(A, b);
        u = sol[:self.num_dof];

        self.rigid_motion_rhs = self.rigid_motion * u;

    def __initialize_history(self):
        self.parameter_history = [];
        self.objective_history = [];
        self.gradient_history = [];
        self.grad_young_history = [];
        self.displacement_history = [];
        self.lagrange_history = [];
        self.displacement_strain_history = [];
        self.lagrange_strain_history = [];
        self.cache = {};
        self.iteration_indices = [];

    def __compute_displacement(self):
        K = self.stiffness;
        R = self.rigid_motion;

        A = scipy.sparse.bmat([[K, R.T], [R, None]]).tocsc();
        b = np.zeros(A.shape[0]);
        b[0:self.num_dof] = self.source_term;
        if hasattr(self, "rigid_motion_rhs"):
            b[self.num_dof:] = self.rigid_motion_rhs;
        else:
            b[self.num_dof:] = R * self.target_displacement;

        sol = scipy.sparse.linalg.spsolve(A, b);
        u = sol[:self.num_dof];

        self.displacement = u;

    def __compute_lagrange_multiplier(self):
        displacement_gap = (self.target_displacement - self.displacement)\
                .reshape((self.mesh.num_vertices, self.mesh.dim), order="C");

        self.assembler.set_material(self.hetero_material);
        K = self.stiffness;
        R = self.rigid_motion;

        self.lagrange_source_term = displacement_gap * self.target_areas[:,np.newaxis];
        self.lagrange_source_term = self.lagrange_source_term.ravel(order="C");

        A = scipy.sparse.bmat([[K, R.T], [R, None]]).tocsc();
        b = np.zeros(A.shape[0]);
        b[0:self.num_dof] = self.lagrange_source_term;

        sol = scipy.sparse.linalg.spsolve(A, b);
        l = sol[:self.num_dof];

        self.lagrange_multiplier = l;

    def __compute_objective_grad(self):
        u_strain = self.displacement_strain * self.displacement;
        l_strain = self.displacement_strain * self.lagrange_multiplier;

        self.assembler.set_material(self.grad_young_material);
        grad_young_elasticity = format(self.assembler.assemble("elasticity_tensor"));
        self.assembler.set_material(self.grad_poisson_material);
        grad_poisson_elasticity = format(self.assembler.assemble("elasticity_tensor"));

        int_u_young_stress = self.element_volumes * grad_young_elasticity * u_strain;
        int_u_poisson_stress = self.element_volumes * grad_poisson_elasticity * u_strain;

        self.grad_young = l_strain * self.strain_doubler * int_u_young_stress;
        self.grad_poisson = l_strain * self.strain_doubler * int_u_poisson_stress;

        self.grad_young = self.grad_young.reshape((self.mesh.num_elements, -1))
        self.grad_poisson = self.grad_poisson.reshape((self.mesh.num_elements, -1))

        self.grad_young = np.sum(self.grad_young, axis=1);
        self.grad_poisson = np.sum(self.grad_poisson, axis=1);

        return self.grad_young, self.grad_poisson;

    def __update_elasticity_gradient(self):
        young = self.mesh.get_attribute(self.young_field_name);
        poisson = self.mesh.get_attribute(self.poisson_field_name);

        grad_young   = [self.C_E(E, v).ravel() for E, v in zip(young, poisson)];
        grad_poisson = [self.C_v(E, v).ravel() for E, v in zip(young, poisson)];

        self.mesh.set_attribute(self.grad_young_field_name, np.hstack(grad_young));
        self.mesh.set_attribute(self.grad_poisson_field_name, np.hstack(grad_poisson));

    def __update_elasticity_model(self):
        self.hetero_material.update();
        self.grad_young_material.update();
        self.grad_poisson_material.update();

        self.assembler.set_material(self.hetero_material);
        self.stiffness = format(self.assembler.assemble("stiffness"));

    def evaluate_regularizer(self, young, poisson):
        young_smoothness = 0.5 * np.dot(young, self.element_connectivity * young);
        poisson_smoothness = 0.5 * np.dot(poisson, self.element_connectivity *
                poisson);
        return young_smoothness + poisson_smoothness;

    def evaluate_regularizer_gradient(self, young, poisson):
        grad_young_smoothness = self.element_connectivity * young;
        grad_poisson_smoothness = self.element_connectivity * poisson;
        return np.hstack((grad_young_smoothness, grad_poisson_smoothness));

    def __evaluate_objective(self):
        displacement_gap = (self.displacement - self.target_displacement)\
                .reshape((self.mesh.num_vertices, self.mesh.dim), order="C");
        #print(self.displacement.reshape((self.mesh.num_vertices,
        #    self.mesh.dim))[self.fixed_idx]);
        obj_value = 0.5 * np.sum(norm(displacement_gap, axis=1)**2 * self.target_areas);
        return obj_value;

    def __evaluate_objective_and_gradient(self, parameters):
        key = self.__parameters_to_key(parameters);
        if key in self.cache:
            idx = self.cache[key];
            return self.objective_history[idx], self.gradient_history[idx];
        else:
            young, poisson = self.__parse_parameters(parameters);
            self.mesh.set_attribute(self.young_field_name, young);
            self.mesh.set_attribute(self.poisson_field_name, poisson);
            self.__update_elasticity_gradient();
            self.__update_elasticity_model();
            self.__compute_displacement();
            self.__compute_lagrange_multiplier();
            self.__compute_objective_grad();

            objective = self.__evaluate_objective()\
                    #+ self.evaluate_regularizer(young, poisson);
            gradient = np.hstack((self.grad_young, self.grad_poisson))\
                    #+ self.evaluate_regularizer_gradient(young, poisson);

            displacement_strain = self.displacement_strain * self.displacement;
            lagrange_strain = self.displacement_strain * self.lagrange_multiplier;

            idx = len(self.cache);
            self.parameter_history.append(np.copy(parameters));
            self.objective_history.append(np.copy(objective));
            self.gradient_history.append(np.copy(gradient));
            self.grad_young_history.append(np.copy(self.grad_young));
            self.displacement_history.append(np.copy(self.displacement));
            self.lagrange_history.append(np.copy(self.lagrange_multiplier));
            self.displacement_strain_history.append(displacement_strain);
            self.lagrange_strain_history.append(lagrange_strain);
            self.cache[key] = idx;

            return objective, gradient;

    def __parse_parameters(self, parameters):
        young = parameters[:self.mesh.num_elements];
        poisson = parameters[self.mesh.num_elements:];
        return young, poisson;

    def __parameters_to_key(self, parameters):
        key = hashlib.sha1(parameters.view(np.uint8)).hexdigest();
        return key;

    @property
    def parameters(self):
        young = self.mesh.get_attribute(self.young_field_name).ravel();
        poisson = self.mesh.get_attribute(self.poisson_field_name).ravel();
        return np.hstack((young, poisson));

    @property
    def bounds(self):
        return [(0.1, 5.0)] * self.mesh.num_elements +\
                [(-0.3, 0.3)] * self.mesh.num_elements;