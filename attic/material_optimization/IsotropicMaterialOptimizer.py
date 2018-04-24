import numpy as np
from numpy.linalg import norm
import scipy.sparse
import scipy.sparse.linalg

from MaterialOptimizer import MaterialOptimizer
import LinearElasticitySettings
import PyAssembler
from timethis import timethis
from MatrixUtils import format

class IsotropicMaterialOptimizer(MaterialOptimizer):
    def __init__(self, mesh):
        super(IsotropicMaterialOptimizer, self).__init__(mesh);
        self.__initialize_material_parameters();
        self.__initialize_elasticity_gradient();
        self.__initialize_materials();
        self.__initialize_assembler();
        self.__initialize_matrices();

    @timethis
    def optimize(self, max_iterations):
        obj_history = [];
        grad_history = [];
        young_history = [];
        for i in range(max_iterations):
            #finite_diff_grad_young = self.__compute_finite_difference_grad();

            self.__compute_displacement();
            self.__compute_lagrange_multiplier();
            grad_young, grad_poisson = self.__compute_objective_grad();
            self.__update_material_parameters(grad_young, grad_poisson);
            self.__update_elasticity_model();

            #grad_difference = norm(finite_diff_grad_young - grad_young);
            obj_val = self.__evaluate_objective();
            grad_norm = norm(grad_young);
            obj_history.append(obj_val);
            grad_history.append(grad_norm);
            young_history.append(np.copy(self.mesh.get_attribute("young")));
            #print("itr: {}  grd: {}  obj: {}  grad_diff: {}".format(
            #    i, grad_norm, obj_val, grad_difference));
            print("itr: {}  grd: {}  obj: {}".format(
                i, grad_norm, obj_val));

            if self.__has_converged(grad_history, obj_history):
                print("Converged in {} iterations.".format(i));
                break;
        self.save_iterations(obj_history, grad_history);
        self.add_young_history_to_mesh(young_history);

    def save_iterations(self, obj_hist, grad_hist):
        assert(len(obj_hist) == len(grad_hist));
        import csv
        log_file = "iteration_history.csv";
        num_data = len(obj_hist);
        with open(log_file, 'w') as fout:
            writer = csv.writer(fout);
            writer.writerow(["iteration", "objective", "gradient"]);
            rows = np.vstack((range(num_data), obj_hist, grad_hist)).T;
            writer.writerows(rows);
        print("iteration history saved to {}".format(log_file));

    def add_young_history_to_mesh(self, young_history):
        for i,young in enumerate(young_history):
            name = "young_{}".format(i);
            self.mesh.add_attribute(name);
            self.mesh.set_attribute(name, young);

    @timethis
    def __initialize_material_parameters(self):
        self.density = 1.0

        self.young_field_name = "young";
        self.mesh.add_attribute(self.young_field_name);
        self.mesh.set_attribute(self.young_field_name,
                5.0 * np.ones(self.mesh.num_elements));

        self.poisson_field_name = "poisson";
        self.mesh.add_attribute(self.poisson_field_name);
        self.mesh.set_attribute(self.poisson_field_name,
                0.0 * np.ones(self.mesh.num_elements));

    @timethis
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
                    self.C_lambda * (E / (1-v**2) + 2*E*v**2/(1-v**2)**2) +\
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

    @timethis
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

    @timethis
    def __initialize_assembler(self):
        self.assembler = PyAssembler.FEAssembler.create(
                self.mesh.raw_mesh, self.hetero_material);

    @timethis
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

    @timethis
    def __update_elasticity_gradient(self):
        young = self.mesh.get_attribute(self.young_field_name);
        poisson = self.mesh.get_attribute(self.poisson_field_name);

        grad_young   = [self.C_E(E, v).ravel() for E, v in zip(young, poisson)];
        grad_poisson = [self.C_v(E, v).ravel() for E, v in zip(young, poisson)];

        self.mesh.set_attribute(self.grad_young_field_name, np.hstack(grad_young));
        self.mesh.set_attribute(self.grad_poisson_field_name, np.hstack(grad_poisson));

    @timethis
    def __update_elasticity_model(self):
        self.hetero_material.update();
        self.grad_young_material.update();
        self.grad_poisson_material.update();

        self.assembler.set_material(self.hetero_material);
        self.stiffness = format(self.assembler.assemble("stiffness"));

    @timethis
    def __update_material_parameters(self, delta_young, delta_poisson):
        young_step_size = 500.0;#abs(0.01 / np.amax(delta_young));
        poisson_step_size = 1.0

        young = self.mesh.get_attribute(self.young_field_name).ravel();
        poisson = self.mesh.get_attribute(self.poisson_field_name).ravel();

        young -= young_step_size * delta_young;
        #poisson -= poisson_step_size * delta_poisson;

        self.mesh.set_attribute(self.young_field_name, young);
        #self.mesh.set_attribute(self.poisson_field_name, poisson);

    @timethis
    def __compute_displacement(self):
        num_dof = self.mesh.dim * self.mesh.num_vertices;

        K = self.stiffness;
        R = self.rigid_motion;

        A = scipy.sparse.bmat([[K, R.T], [R, None]]).tocsc();
        b = np.zeros(A.shape[0]);
        b[0:num_dof] = self.source_term;
        b[num_dof:] = R * self.target_displacement;

        sol = scipy.sparse.linalg.spsolve(A, b);
        u = sol[:num_dof];

        self.displacement = u;

    @timethis
    def __compute_lagrange_multiplier(self):
        num_dof = self.mesh.dim * self.mesh.num_vertices;
        displacement_gap = (self.target_displacement - self.displacement)\
                .reshape((self.mesh.num_vertices, self.mesh.dim), order="C");

        self.assembler.set_material(self.hetero_material);
        K = self.stiffness;
        R = self.rigid_motion;

        self.lagrange_source_term = displacement_gap * self.target_areas[:,np.newaxis];
        self.lagrange_source_term = self.lagrange_source_term.ravel(order="C");

        A = scipy.sparse.bmat([[K, R.T], [R, None]]).tocsc();
        b = np.zeros(A.shape[0]);
        b[0:num_dof] = self.lagrange_source_term;

        sol = scipy.sparse.linalg.spsolve(A, b);
        l = sol[:num_dof];

        self.lagrange_multiplier = l;

    @timethis
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

    @timethis
    def __compute_finite_difference_grad(self):
        def get_objective_value(young_array):
            self.mesh.set_attribute("young", young_array);
            self.hetero_material.update();
            self.assembler.set_material(self.hetero_material);
            self.stiffness = format(self.assembler.assemble("stiffness"));

            self.__compute_displacement();
            return self.__evaluate_objective();

        young = np.copy(self.mesh.get_attribute("young"));
        step_size = 1e-3;
        ori_obj = get_objective_value(young);

        young_grad = np.zeros(self.mesh.num_elements);
        for i in range(self.mesh.num_elements):
            cur_young = np.copy(young);
            cur_young[i] += step_size;
            cur_obj = get_objective_value(cur_young);
            young_grad[i] = (cur_obj - ori_obj) / step_size;

        self.mesh.set_attribute("young", young);
        self.hetero_material.update();
        self.stiffness = format(self.assembler.assemble("stiffness"));
        return young_grad;

    @timethis
    def __evaluate_objective(self):
        displacement_gap = (self.displacement - self.target_displacement)\
                .reshape((self.mesh.num_vertices, self.mesh.dim), order="C");
        obj_value = 0.5 * np.sum(norm(displacement_gap, axis=1)**2 * self.target_areas);
        return obj_value;

    def __has_converged(self, obj_hist, grad_hist):
        return grad_hist[-1] < grad_hist[0] * 1e-2;

