import LinearElasticitySettings

import numpy as np
from numpy.linalg import norm

from Assembler import Assembler
from Dual import Dual
from Primal import Primal
from MaterialState import MaterialState

class MaterialOptimizer(object):
    def __init__(self, mesh, dirichlet_bc, neumann_bc):
        self.mesh = mesh;
        self.dirichlet_bc = dirichlet_bc;
        self.neumann_bc = neumann_bc;

        self.material = MaterialState(mesh);
        self.assembler = Assembler(mesh, self.material);

        self.primal = Primal(mesh, self.assembler);
        self.primal.add_dirichlet_bc(dirichlet_bc);
        self.primal.add_neumann_bc(neumann_bc);

        self.dual = Dual(mesh, self.assembler);
        self.dual.add_neumann_bc(neumann_bc);
        self.init_log();

        self.objective_tol = 1e-3;
        self.sufficient_decrease_tol = 0.05;

    def optimize(self, max_iter = 10):
        self.iterations = 0;
        for i in range(max_iter):
            self.assembler.update();
            self.primal.solve();
            self.dual.set_rigid_motion_as(self.primal.displacement);
            self.dual.solve();
            self.fit_material();
            self.iterations += 1;

            self.evaluate_objective();
            self.log_progress();

            print("{:3}: obj={}".format(i, self.objective_value));
            if self.has_converged():
                break;

    def fit_material(self):
        num_elements = self.mesh.num_elements;
        strain = self.primal.strain.reshape((num_elements, -1), order="C");
        stress = self.dual.stress.reshape((num_elements, -1), order="C");

        young = np.zeros(self.mesh.num_elements);
        poisson = np.zeros(self.mesh.num_elements);
        for i, epsilon, sigma in zip(range(num_elements), strain, stress):
            # So far, just fit Young's modulus.
            # Assume poisson is zero.
            young[i] = norm(sigma) / norm(epsilon);

        self.material.update(young, poisson);

    def evaluate_objective(self):
        u_p = self.primal.displacement;
        u_d = self.dual.displacement;
        u_gap = (u_p - u_d).reshape((self.mesh.num_nodes, -1), order="C");
        self.objective_value = np.dot(
                np.square(norm(u_gap[self.dirichlet_bc[0]], axis=1)),
                self.dirichlet_bc[2])

    def init_log(self):
        self.primal_displacement = [];
        self.dual_displacement = [];
        self.young = [];
        self.objective_history = [];

    def log_progress(self):
        self.primal_displacement.append(np.copy(self.primal.displacement));
        self.dual_displacement.append(np.copy(self.dual.displacement));
        self.young.append(np.copy(self.mesh.get_attribute(self.material.young_attr_name)));
        self.objective_history.append(self.objective_value);

    def has_converged(self):
        if self.objective_value < self.objective_tol:
            print("Converged because objective {} < {}".format(
                self.objective_value, self.objective_tol));
            return True;
        if len(self.objective_history) > 1:
            prev_obj = self.objective_history[-2];
            curr_obj = self.objective_history[-1];
            rel_decrease = (prev_obj - curr_obj) / prev_obj;
            if rel_decrease < self.sufficient_decrease_tol:
                print("Converged because relative objective decrease is {} < {}"\
                        .format(rel_decrease, self.sufficient_decrease_tol));
                return True;
        return False;

