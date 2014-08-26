import json
import os.path

import numpy as np

class MaterialParameter(object):
    def __init__(self, dim, param_file):
        self.dim = dim;
        if dim == 2:
            self.__extract_2D_param(param_file);
        elif dim == 3:
            self.__extract_3D_param(param_file);

    def __extract_2D_param(self, param_file):
        with open(param_file, 'r') as fin:
            self.config = json.load(fin);

        self.young_x = self.config["youngs_modulus"][0];
        self.young_y = self.config["youngs_modulus"][1];
        self.shear_xy = self.config["shear_modulus"][0];
        self.poisson_xy = self.config["poisson_ratio"][0];
        self.poisson_yx = self.config["poisson_ratio"][1];

        self.names = [
                "young_x",
                "young_y",
                "shear_xy",
                "poisson_xy",
                "poisson_yx" ];

    def __extract_3D_param(self, param_file):
        with open(param_file, 'r') as fin:
            self.config = json.load(fin);

        self.young_x = self.config["youngs_modulus"][0];
        self.young_y = self.config["youngs_modulus"][1];
        self.young_z = self.config["youngs_modulus"][2];
        self.shear_yz = self.config["shear_modulus"][0];
        self.shear_zx = self.config["shear_modulus"][1];
        self.shear_xy = self.config["shear_modulus"][2];
        self.poisson_yz = self.config["poisson_ratio"][0];
        self.poisson_zy = self.config["poisson_ratio"][1];
        self.poisson_zx = self.config["poisson_ratio"][2];
        self.poisson_xz = self.config["poisson_ratio"][3];
        self.poisson_xy = self.config["poisson_ratio"][4];
        self.poisson_yx = self.config["poisson_ratio"][5];

        self.names = [
                "young_x",
                "young_y",
                "young_z",
                "shear_yz",
                "shear_zx",
                "shear_xy",
                "poisson_yz",
                "poisson_zy",
                "poisson_zx",
                "poisson_xz",
                "poisson_xy",
                "poisson_yx",
                ];

    @property
    def names(self):
        return self.__names;

    @names.setter
    def names(self, value):
        self.__names = value;

    @property
    def values(self):
        return [getattr(self, name) for name in self.names];

