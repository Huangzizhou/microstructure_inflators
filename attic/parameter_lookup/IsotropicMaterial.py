from Material import Material

import numpy as np
from numpy.linalg import inv

class IsotropicMaterial(Material):
    def __init__(self, dim, young, poisson):
       self.dim = dim; 
       self.young = young;
       self.poisson = poisson;
       if self.dim == 2:
           self.__generate_tensor_2D();
       elif self.dim == 3:
           self.__generate_tensor_3D();
       else:
            raise NotImplementedError("Dim {} is not supported.".format(dim));

    def __generate_tensor_2D(self):
        """
        Formula from
        http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm
        """
        E = self.young;
        nu = self.poisson;
        self.compliance_tensor = np.array([
            [1.0/E, -nu/E, 0.0],
            [-nu/E, 1.0/E, 0.0],
            [0.0, 0.0, (1.0+nu)/E]
            ]);
        self.elasticity_tensor = inv(self.compliance_tensor);

    def __generate_tensor_3D(self):
        """
        Formula from
        http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm
        """
        E = self.young;
        nu = self.poisson;
        self.compliance_tensor = np.array([
            [1.0/E, -nu/E, -nu/E, 0.0, 0.0, 0.0],
            [-nu/E, 1.0/E, -nu/E, 0.0, 0.0, 0.0],
            [-nu/E, -nu/E, 1.0/E, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 2*(1.0+nu)/E, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 2*(1.0+nu)/E, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 2*(1.0+nu)/E] ]);
        self.elasticity_tensor = inv(self.compliance_tensor);
