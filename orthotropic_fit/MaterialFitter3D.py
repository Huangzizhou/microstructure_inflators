import numpy as np
from numpy.linalg import lstsq, inv, norm

import LinearElasticitySettings
from MaterialFitter import MaterialFitter
import PredefinedBoundaryConditions as PredefinedBC
from ElasticityUtils import displacement_to_strain
from timethis import timethis

class MaterialFitter3D(MaterialFitter):
    def __init__(self, mesh, material_file=None):
        super(MaterialFitter3D, self).__init__(mesh, material_file);

    @property
    def bc_configs(self):
        bbox_min, bbox_max = self.mesh.bbox;
        eps = min(2e-1, norm(bbox_max - bbox_min) * 0.01);
        return [PredefinedBC.compress_x(bbox_min, bbox_max, eps),
                PredefinedBC.compress_y(bbox_min, bbox_max, eps),
                PredefinedBC.compress_z(bbox_min, bbox_max, eps),
                PredefinedBC.compress_xy(bbox_min, bbox_max, eps),
                PredefinedBC.compress_yz(bbox_min, bbox_max, eps),
                PredefinedBC.compress_zx(bbox_min, bbox_max, eps),
                PredefinedBC.shear_xy(bbox_min, bbox_max, eps),
                PredefinedBC.shear_yz(bbox_min, bbox_max, eps),
                PredefinedBC.shear_zx(bbox_min, bbox_max, eps) ];

    @timethis
    def _fit_material_parameters(self):
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

                coeff_44 = np.multiply(strain_1_xy, strain_2_xy).dot(voxel_volumes) * 4;
                coeff_55 = np.multiply(strain_1_xz, strain_2_xz).dot(voxel_volumes) * 4;
                coeff_66 = np.multiply(strain_1_yz, strain_2_yz).dot(voxel_volumes) * 4;

                coeff_12 = (np.multiply(strain_1_xx, strain_2_yy) +\
                        np.multiply(strain_1_yy, strain_2_xx)).dot(voxel_volumes)
                coeff_13 = (np.multiply(strain_1_xx, strain_2_zz) +\
                        np.multiply(strain_1_zz, strain_2_xx)).dot(voxel_volumes)
                coeff_23 = (np.multiply(strain_1_yy, strain_2_zz) +\
                        np.multiply(strain_1_zz, strain_2_yy)).dot(voxel_volumes)

                coeff_14 = (np.multiply(strain_1_xx, strain_2_xy) +\
                        np.multiply(strain_1_xy, strain_2_xx)).dot(voxel_volumes) * 2;
                coeff_15 = (np.multiply(strain_1_xx, strain_2_xz) +\
                        np.multiply(strain_1_xz, strain_2_xx)).dot(voxel_volumes) * 2;
                coeff_16 = (np.multiply(strain_1_xx, strain_2_yz) +\
                        np.multiply(strain_1_yz, strain_2_xx)).dot(voxel_volumes) * 2;

                coeff_24 = (np.multiply(strain_1_yy, strain_2_xy) +\
                        np.multiply(strain_1_xy, strain_2_yy)).dot(voxel_volumes) * 2;
                coeff_25 = (np.multiply(strain_1_yy, strain_2_xz) +\
                        np.multiply(strain_1_xz, strain_2_yy)).dot(voxel_volumes) * 2;
                coeff_26 = (np.multiply(strain_1_yy, strain_2_yz) +\
                        np.multiply(strain_1_yz, strain_2_yy)).dot(voxel_volumes) * 2;

                coeff_34 = (np.multiply(strain_1_zz, strain_2_xy) +\
                        np.multiply(strain_1_xy, strain_2_zz)).dot(voxel_volumes) * 2;
                coeff_35 = (np.multiply(strain_1_zz, strain_2_xz) +\
                        np.multiply(strain_1_xz, strain_2_zz)).dot(voxel_volumes) * 2;
                coeff_36 = (np.multiply(strain_1_zz, strain_2_yz) +\
                        np.multiply(strain_1_yz, strain_2_zz)).dot(voxel_volumes) * 2;

                coeff_45 = (np.multiply(strain_1_xz, strain_2_xy) +\
                        np.multiply(strain_1_xy, strain_2_xz)).dot(voxel_volumes) * 4;
                coeff_46 = (np.multiply(strain_1_yz, strain_2_xy) +\
                        np.multiply(strain_1_xy, strain_2_yz)).dot(voxel_volumes) * 4;
                coeff_56 = (np.multiply(strain_1_yz, strain_2_xz) +\
                        np.multiply(strain_1_xz, strain_2_yz)).dot(voxel_volumes) * 4;

                coeff.append([
                    coeff_11, coeff_22, coeff_33,
                    coeff_44, coeff_55, coeff_66,
                    coeff_12, coeff_13, coeff_23,
                    coeff_14, coeff_15, coeff_16,
                    coeff_24, coeff_25, coeff_26,
                    coeff_34, coeff_35, coeff_36,
                    coeff_45, coeff_46, coeff_56]);
                rhs.append(energy);

        parameter, residual, rank, singular_vals =\
                lstsq(coeff, rhs);
        self.residual_error = norm((rhs - np.dot(coeff, parameter)) / rhs);
        C_idx = np.array([[ 0, 6, 7, 9,10,11],
                          [ 6, 1, 8,12,13,14],
                          [ 7, 8, 2,15,16,17],
                          [ 9,12,15, 3,18,19],
                          [10,13,16,18, 4,20],
                          [11,14,17,19,20, 5] ]);
        C = parameter[C_idx];
        self.elasticity_tensor = C;

        S = inv(C);
        self.compliance_tensor = S;
        self.condition_num = np.max(singular_vals) / np.min(singular_vals);

