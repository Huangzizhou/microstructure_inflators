import unittest
import numpy as np
from numpy.linalg import norm

import LinearElasticitySettings
from mesh_io import load_mesh, save_mesh
from BoundaryConditionExtractor import BoundaryConditionExtractor
from IsotropicMatOptSetting import IsotropicMatOptSetting

class IsotropicMatOptSettingTest(unittest.TestCase):
    def setUp(self):
        self.mesh = load_mesh("example/square_4x4.obj");
        bc = self.extract_bc("example/tilt.bc");
        self.opt_setting = IsotropicMatOptSetting(self.mesh, bc,
                np.ones(self.mesh.num_elements),
                np.zeros(self.mesh.num_elements));

    def extract_bc(self, bd_file):
        bc_extractor = BoundaryConditionExtractor(self.mesh);
        bc_extractor.extract_from_file(bd_file);
        return bc_extractor;

    def finite_diff_objective_grad(self, parameters):
        epsilon = 1e-6;
        grad = np.zeros(len(parameters));
        ori_obj, ori_grad = self.opt_setting.evaluate(parameters);
        for i in range(len(parameters)):
            probe = np.copy(parameters);
            probe[i] += epsilon;
            cur_obj, cur_grad = self.opt_setting.evaluate(probe);
            grad[i] = (cur_obj - ori_obj) / epsilon;

        return grad;

    def finite_diff_regularizer_grad(self, parameters):
        epsilon = 1e-6;
        young = parameters[:self.mesh.num_elements];
        poisson = parameters[self.mesh.num_elements:];

        value = self.opt_setting.evaluate_regularizer(young, poisson);

        grad = np.zeros(len(parameters));
        for i in range(self.mesh.num_elements):
            mod_young = np.copy(young);
            mod_young[i] += epsilon;
            mod_poisson = np.copy(poisson);
            mod_poisson[i] += epsilon;

            value_mod_young = self.opt_setting.evaluate_regularizer(
                    mod_young, poisson);
            value_mod_poisson = self.opt_setting.evaluate_regularizer(
                    young, mod_poisson);
            grad[i] = (value_mod_young - value) / epsilon;
            grad[i+self.mesh.num_elements] = (value_mod_poisson - value) / epsilon;
        return grad;

    def finite_difference_elasticity_matrix_gradient(self, parameters):
        C = lambda E,v : E/(1-v**2) * np.array(
                [ [1.0, v, 0], [v, 1.0, 0.0], [0.0, 0.0, 0.5*(1.0-v)] ]);
        young = parameters[:self.mesh.num_elements];
        poisson = parameters[self.mesh.num_elements:];

        epsilon = 1e-6;
        C_grad_E = [];
        C_grad_v = [];
        for E,v in zip(young, poisson):
            ori_C = C(E, v);
            mod_E_C = C(E+epsilon, v);
            mod_v_C = C(E, v+epsilon);
            dC_dE = (mod_E_C - ori_C) / epsilon;
            dC_dv = (mod_v_C - ori_C) / epsilon;
            C_grad_E.append(dC_dE);
            C_grad_v.append(dC_dv);
        return C_grad_E, C_grad_v;

    def add_attribute(self, name, val):
        self.mesh.add_attribute(name);
        self.mesh.set_attribute(name, val);

    def get_parameters(self):
        num_parameters = self.mesh.num_elements * 2;
        parameters = np.ones(num_parameters);
        for i in range(self.mesh.num_elements):
            parameters[i] = 5.0 + 0.1 * i;
            parameters[i + self.mesh.num_elements] = 0.3 + 0.01 * i;
        return parameters;

    def test_initialization(self):
        self.assertTrue(self.mesh.has_attribute(
            self.opt_setting.young_field_name));
        self.assertTrue(self.mesh.has_attribute(
            self.opt_setting.poisson_field_name));
        self.assertTrue(self.mesh.has_attribute(
            self.opt_setting.grad_young_field_name));
        self.assertTrue(self.mesh.has_attribute(
            self.opt_setting.grad_poisson_field_name));

    def test_initial_parameter(self):
        young_value = 1.5;
        poisson_value = 0.3;
        num_parameters = self.mesh.num_elements * 2;
        parameters = np.ones(num_parameters) * young_value;
        parameters[self.mesh.num_elements:] = poisson_value;

        obj, grad = self.opt_setting.evaluate(parameters);

        young = self.mesh.get_attribute(self.opt_setting.young_field_name);
        poisson = self.mesh.get_attribute(self.opt_setting.poisson_field_name);

        self.assertEqual(young_value, np.amax(young));
        self.assertEqual(young_value, np.amin(young));
        self.assertEqual(poisson_value, np.amax(poisson));
        self.assertEqual(poisson_value, np.amin(poisson));

    def test_finite_diff_elasticity_matrix_grad(self):
        dim = self.mesh.dim;
        tensor_size = dim * (dim+1) / 2;
        parameters = self.get_parameters();
        obj, grad = self.opt_setting.evaluate(parameters);

        grad_C_E = self.mesh.get_attribute(self.opt_setting.grad_young_field_name);
        grad_C_v = self.mesh.get_attribute(self.opt_setting.grad_poisson_field_name);
        grad_C_E = grad_C_E.reshape((self.mesh.num_elements, -1), order="C");
        grad_C_v = grad_C_v.reshape((self.mesh.num_elements, -1), order="C");

        fd_grad_C_E, fd_grad_C_v = \
                self.finite_difference_elasticity_matrix_gradient(parameters);

        for i in range(self.mesh.num_elements):
            grad_C_E_i = grad_C_E[i].reshape((tensor_size, tensor_size));
            grad_C_v_i = grad_C_v[i].reshape((tensor_size, tensor_size));
            fd_grad_C_E_i = fd_grad_C_E[i];
            fd_grad_C_v_i = fd_grad_C_v[i];

            self.assertAlmostEqual(0.0, norm(grad_C_E_i - fd_grad_C_E_i), 3);
            self.assertAlmostEqual(0.0, norm(grad_C_v_i - fd_grad_C_v_i), 3);

    def test_regularizer_evaluation(self):
        parameters = self.get_parameters();
        regularizer_grad = self.opt_setting.evaluate_regularizer_gradient(
                parameters[:self.mesh.num_elements],
                parameters[self.mesh.num_elements:]);
        regularizer_finite_diff_grad = self.finite_diff_regularizer_grad(parameters);
        diff_grad = regularizer_grad - regularizer_finite_diff_grad;
        self.assertAlmostEqual(0.0, norm(diff_grad), 3);

    #@unittest.skip("debugging")
    def test_obj_evaluation(self):
        parameters = self.get_parameters();
        obj, grad = self.opt_setting.evaluate(parameters);

        self.add_attribute("displacement",
                self.opt_setting.displacement);
        self.add_attribute("target_displacement",
                self.opt_setting.target_displacement);
        self.add_attribute("multiplier",
                self.opt_setting.lagrange_multiplier);
        save_mesh("tmp.msh", self.mesh,
                self.opt_setting.young_field_name,
                self.opt_setting.poisson_field_name,
                "displacement",
                "target_displacement",
                "multiplier"
                );

        diff_grad = self.finite_diff_objective_grad(parameters);

        self.assertAlmostEqual(0.0, norm((grad - diff_grad)), 3);

