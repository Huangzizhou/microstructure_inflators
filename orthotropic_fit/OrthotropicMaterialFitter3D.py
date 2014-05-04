import numpy as np
from numpy.linalg import lstsq, inv, norm

import LinearElasticitySettings

from MaterialFitter3D import MaterialFitter3D

from ElasticityUtils import displacement_to_strain
from timethis import timethis

class OrthotropicMaterialFitter3D(MaterialFitter3D):
    def __init__(self, mesh, material_file=None):
        super(OrthotropicMaterialFitter3D, self).__init__(mesh, material_file);

    @timethis
    def _fit_material_parameters(self):
        return self.__fit_orthotropic_parameters_3D();

    @timethis
    def __fit_orthotropic_parameters_3D(self):
        tensor_size = 6;

        voxel_volumes = self.coarse_mesh.get_attribute("voxel_volume").ravel();
        K = self.assembler.stiffness_matrix;

        coeff = [];
        rhs = [];
        for u1,coarse_u1 in zip(self.displacements, self.coarse_displacements):
            strain_1 = displacement_to_strain(self.coarse_assembler, coarse_u1);
            strain_1 = strain_1.reshape((-1, tensor_size), order="C");

            strain_1_xx = strain_1[:,0].ravel();
            strain_1_yy = strain_1[:,1].ravel();
            strain_1_zz = strain_1[:,2].ravel();
            strain_1_xy = strain_1[:,3].ravel();
            strain_1_xz = strain_1[:,4].ravel();
            strain_1_yz = strain_1[:,5].ravel();
            for u2,coarse_u2 in zip(self.displacements, self.coarse_displacements):
                energy = np.dot(u1, K*u2);

                strain_2 = displacement_to_strain(self.coarse_assembler, coarse_u2);
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
                coeff_44 = np.multiply(strain_1_yz, strain_2_yz).dot(voxel_volumes) * 4;
                coeff_55 = np.multiply(strain_1_xz, strain_2_xz).dot(voxel_volumes) * 4;
                coeff_66 = np.multiply(strain_1_xy, strain_2_xy).dot(voxel_volumes) * 4;
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
        self.elasticity_tensor = C;

        S = inv(C);
        self.compliance_tensor = S;
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


    @property
    def youngs_modulus(self):
        dim = self.mesh.dim;
        young = 1.0 / self.orthotropic_parameter[:dim];
        if np.any(young <= 0):
            import warnings
            warnings.warn("Negative Young's modulus: {}".format(young));
        return abs(young);

    @property
    def poisson_ratio(self):
        dim = self.mesh.dim;
        young = self.youngs_modulus;
        # Order is important here!
        # It is consistent with the ordering of input for Orthotropic
        # material of PyAssembler.
        return np.array( [
            -self.orthotropic_parameter[8] * young[1], # v_yz
            -self.orthotropic_parameter[8] * young[2], # v_zy

            -self.orthotropic_parameter[7] * young[2], # v_zx
            -self.orthotropic_parameter[7] * young[0], # v_xz

            -self.orthotropic_parameter[6] * young[0], # v_xy
            -self.orthotropic_parameter[6] * young[1], # v_yx
            ]);

    @property
    def shear_modulus(self):
        dim = self.mesh.dim;
        # Order: [G_yz, G_zx, G_xy]
        shear = 1.0 / self.orthotropic_parameter[3:6];

        if np.any(shear <= 0):
            import warnings
            warnings.warn("Negative shear modulus: {}".format(shear));
        return abs(shear);

