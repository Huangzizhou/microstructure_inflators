import numpy as np
from numpy.linalg import lstsq, inv, norm

import LinearElasticitySettings

from MaterialFitter2D import MaterialFitter2D

from ElasticityUtils import displacement_to_strain
from timethis import timethis

class OrthotropicMaterialFitter2D(MaterialFitter2D):
    def __init__(self, mesh, material_file=None):
        super(OrthotropicMaterialFitter2D, self).__init__(mesh, material_file);

    @timethis
    def _fit_material_parameters(self):
        return self.__fit_orthotropic_parameters_2D();

    @timethis
    def __fit_orthotropic_parameters_2D(self):
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

                coeff.append([coeff_11, coeff_22, coeff_33, coeff_12]);
                rhs.append(energy);

        parameter, residual, rank, singular_vals =\
                lstsq(coeff, rhs);
        self.residual_error = norm((rhs - np.dot(coeff, parameter)) / rhs);
        C = np.array([
            [parameter[0], parameter[3],          0.0],
            [parameter[3], parameter[1],          0.0],
            [         0.0,          0.0, parameter[2]] ]);
        self.elasticity_tensor = C;
        S = inv(C);
        self.compliance_tensor = S;
        parameter = [S[0,0], S[1,1], S[2,2], S[0,1]];
        self.orthotropic_parameter = np.array(parameter);
        self.condition_num = np.max(singular_vals) / np.min(singular_vals);


    @property
    def youngs_modulus(self):
        dim = self.mesh.dim;
        young = 1.0 / self.orthotropic_parameter[:dim];
        if np.any(young <= 0.0):
            import warnings
            warnings.warn("Negative Young's modulus: {}".format(young));
        return abs(young);

    @property
    def poisson_ratio(self):
        dim = self.mesh.dim;
        young = self.youngs_modulus;
        return np.array([
            -self.orthotropic_parameter[3] * young[0], # v_xy
            -self.orthotropic_parameter[3] * young[1]  # v_yx
            ]);

    @property
    def shear_modulus(self):
        dim = self.mesh.dim;
        # Order: [G_xy]
        shear = 1.0 / self.orthotropic_parameter[2:3];

        if np.any(shear <= 0):
            import warnings
            warnings.warn("Negative shear modulus: {}".format(shear));
        return abs(shear);

