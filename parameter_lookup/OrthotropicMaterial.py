from Material import Material

import numpy as np
from numpy.linalg import inv

class OrthotropicMaterial(Material):
    def __init__(self, young, poisson, shear):
        if (len(young) == 2):
            assert(len(poisson) == 2);
            assert(len(shear) == 1);
            self.dim = 2;
            self.young = young;
            self.poisson = poisson;
            self.shear = shear;
            self.__generate_tensors_2D();
        elif (len(young) == 3):
            assert(len(poisson) == 6);
            assert(len(shear) == 3);
            self.dim = 3;
            self.young = young;
            self.poisson = poisson;
            self.shear = shear;
            self.__generate_tensors_3D();
        else:
            raise NotImplementedError("Dim {} is not supported.".format(
                len(young)));

    def __generate_tensors_2D(self):
        """
        Formula from
        http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm
        """
        E_x,E_y = self.young;
        nu_xy,nu_yx = self.poisson;
        G_xy = self.shear[0];

        self.compliance_tensor = np.array([
            [1.0/E_x, -nu_yx/E_y, 0.0,],
            [-nu_xy/E_x, 1.0/E_y, 0.0,],
            [0.0, 0.0, 1.0/G_xy] ]);
        self.elasticity_tensor = inv(self.compliance_tensor);

    def __generate_tensors_3D(self):
        """
        Formula from
        http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm
        """
        E_x,E_y,E_z = self.young;
        nu_yz,nu_zy,nu_zx,nu_xz,nu_xy,nu_yz = self.poisson;
        G_yz,G_zx,G_xy = self.shear;

        self.compliance_tensor = np.array([
            [1.0/E_x, -nu_yx/E_y, -nu_zx/E_z, 0.0, 0.0, 0.0,],
            [-nu_xy/E_x, 1.0/E_y, -nu_zy/E_z, 0.0, 0.0, 0.0,],
            [-nu_xz/E_x, -nu_yz/E_y, 1.0/E_z, 0.0, 0.0, 0.0,],
            [0.0, 0.0, 0.0, 1.0/G_yz, 0.0, 0.0,],
            [0.0, 0.0, 0.0, 0.0, 1.0/G_zx, 0.0,],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0/G_xy,] ]);
        self.elasticity_tensor = inv(self.compliance_tensor);
