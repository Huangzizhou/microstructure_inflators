from numpy.linalg import norm

import LinearElasticitySettings
from MaterialFitter import MaterialFitter
import PredefinedBoundaryConditions as PredefinedBC

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
                PredefinedBC.compress_zx(bbox_min, bbox_max, eps) ];

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
                coeff_44 = np.multiply(strain_1_yz, strain_2_yz).dot(voxel_volumes) * 2;
                coeff_55 = np.multiply(strain_1_xz, strain_2_xz).dot(voxel_volumes) * 2;
                coeff_66 = np.multiply(strain_1_xy, strain_2_xy).dot(voxel_volumes) * 2;
                coeff_12 = (np.multiply(strain_1_xx, strain_2_yy) +\
                        np.multiply(strain_1_yy, strain_2_xx)).dot(voxel_volumes)
                coeff_13 = (np.multiply(strain_1_xx, strain_2_zz) +\
                        np.multiply(strain_1_zz, strain_2_xx)).dot(voxel_volumes)
                coeff_23 = (np.multiply(strain_1_yy, strain_2_zz) +\
                        np.multiply(strain_1_zz, strain_2_yy)).dot(voxel_volumes)
                coeff_14 = (np.multiply(strain_1_xx, strain_2_yz) +\
                        np.multiply(strain_1_yz, strain_2_xx)).dot(voxel_volumes) * 2;

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

