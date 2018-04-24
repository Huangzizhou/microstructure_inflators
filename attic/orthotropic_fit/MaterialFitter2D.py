import numpy as np
from numpy.linalg import lstsq, inv, norm

import LinearElasticitySettings
from MaterialFitter import MaterialFitter
import PredefinedBoundaryConditions as PredefinedBC

from ElasticityUtils import displacement_to_strain
from timethis import timethis

class MaterialFitter2D(MaterialFitter):
    def __init__(self, mesh, material_file=None):
        super(MaterialFitter2D, self).__init__(mesh, material_file);

    @property
    def bc_configs(self):
        eps = 1e-3;
        bbox_min, bbox_max = self.mesh.bbox;
        return [PredefinedBC.compress_x(bbox_min, bbox_max, eps),
                PredefinedBC.compress_y(bbox_min, bbox_max, eps),
                PredefinedBC.shear_xy(bbox_min, bbox_max, eps)];

    @timethis
    def _fit_material_parameters(self):
        tensor_size = 3;

        face_areas = self.coarse_mesh.get_attribute("face_area").ravel();
        K = self.assembler.stiffness_matrix;

        coeff = [];
        rhs = [];
        for u1,coarse_u1 in zip(self.displacements, self.coarse_displacements):
            strain_1 = displacement_to_strain(self.coarse_assembler, coarse_u1);
            strain_1 = strain_1.reshape((-1, tensor_size), order="C");

            strain_1_xx = strain_1[:,0].ravel();
            strain_1_yy = strain_1[:,1].ravel();
            strain_1_xy = strain_1[:,2].ravel();

            for u2,coarse_u2 in zip(self.displacements, self.coarse_displacements):
                energy = np.dot(u1, K*u2);

                strain_2 = displacement_to_strain(self.coarse_assembler, coarse_u2);
                strain_2 = strain_2.reshape((-1, tensor_size), order="C");

                strain_2_xx = strain_2[:,0].ravel();
                strain_2_yy = strain_2[:,1].ravel();
                strain_2_xy = strain_2[:,2].ravel();

                coeff_11 = np.multiply(strain_1_xx, strain_2_xx).dot(face_areas);
                coeff_22 = np.multiply(strain_1_yy, strain_2_yy).dot(face_areas);
                coeff_33 = np.multiply(strain_1_xy, strain_2_xy).dot(face_areas) * 4;
                coeff_12 = (np.multiply(strain_1_xx, strain_2_yy) + \
                        np.multiply(strain_1_yy, strain_2_xx)).dot(face_areas);
                coeff_13 = (np.multiply(strain_1_xx, strain_2_xy) +\
                        np.multiply(strain_1_xy, strain_2_xx)).dot(face_areas) * 2;
                coeff_23 = (np.multiply(strain_1_yy, strain_2_xy) +\
                        np.multiply(strain_1_xy, strain_2_yy)).dot(face_areas) * 2;

                coeff.append([coeff_11, coeff_22, coeff_33,
                    coeff_12, coeff_13, coeff_23]);
                rhs.append(energy);

        parameter, residual, rank, singular_vals =\
                lstsq(coeff, rhs);
        self.residual_error = norm((rhs - np.dot(coeff, parameter)) / rhs);
        C_idx = np.array([
            [0, 3, 4],
            [3, 1, 5],
            [4, 5, 2] ]);
        C = parameter[C_idx];
        self.elasticity_tensor = C;
        S = inv(C);
        self.compliance_tensor = S;

        self.condition_num = np.max(singular_vals) / np.min(singular_vals);

