import json
import os.path

import numpy as np
from OrthotropicMaterial import OrthotropicMaterial

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
        self.elasticity_mode_0 = self.config["elasticity_modes"][0];
        self.elasticity_mode_1 = self.config["elasticity_modes"][1];
        self.elasticity_mode_2 = self.config["elasticity_modes"][2];

        if self.poisson_xy > 1.0 or self.poisson_xy < -1.0 or\
                self.poisson_yx > 1.0 or self.poisson_yx < -1.0:
            print("In {}, (v_xy={} v_yx={}) outside of [-1, 1]".format(
                param_file, self.poisson_xy, self.poisson_yx));

        self.names = [
                "young_x",
                "young_y",
                "shear_xy",
                "poisson_xy",
                "poisson_yx",
                "elasticity_mode_0",
                "elasticity_mode_1",
                "elasticity_mode_0",
                ];
        self.material = OrthotropicMaterial(
                [self.young_x, self.young_y],
                [self.poisson_xy, self.poisson_yx],
                [self.shear_xy]);

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
        self.elasticity_mode_0 = self.config["elasticity_modes"][0];
        self.elasticity_mode_1 = self.config["elasticity_modes"][1];
        self.elasticity_mode_2 = self.config["elasticity_modes"][2];
        self.elasticity_mode_3 = self.config["elasticity_modes"][3];
        self.elasticity_mode_4 = self.config["elasticity_modes"][4];
        self.elasticity_mode_5 = self.config["elasticity_modes"][5];

        #if self.poisson_xy > 0.5 or self.poisson_xy < -1.0 or\
        #        self.poisson_yx > 1.0 or self.poisson_yx < -1.0:
        #    print("In {}, (v_xy={} v_yx={}) outside of [-1, 0.5]".format(
        #        param_file, self.poisson_xy, self.poisson_yx));

        #if self.poisson_yz > 0.5 or self.poisson_yz < -1.0 or\
        #        self.poisson_zy > 1.0 or self.poisson_zy < -1.0:
        #    print("In {}, (v_yz={} v_zy={}) outside of [-1, 0.5]".format(
        #        param_file, self.poisson_yz, self.poisson_zy));

        #if self.poisson_zx > 0.5 or self.poisson_zx < -1.0 or\
        #        self.poisson_xz > 1.0 or self.poisson_xz < -1.0:
        #    print("In {}, (v_zx={} v_xz={}) outside of [-1, 0.5]".format(
        #        param_file, self.poisson_zx, self.poisson_xz));

        if self.elasticity_mode_0 / self.elasticity_mode_1 > 100.0:
            print("Pentamode material: {}".format(param_file));

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
                "elasticity_mode_0",
                "elasticity_mode_1",
                "elasticity_mode_2",
                "elasticity_mode_3",
                "elasticity_mode_4",
                "elasticity_mode_5",
                ];

        self.material = OrthotropicMaterial(
                [self.young_x, self.young_y, self.young_z],
                [self.poisson_yz, self.poisson_zy,
                    self.poisson_zx, self.poisson_xz,
                    self.poisson_xy, self.poisson_yx],
                [self.shear_yz, self.shear_zx, self.shear_xy]);

    @property
    def names(self):
        return self.__names;

    @names.setter
    def names(self, value):
        self.__names = value;

    @property
    def values(self):
        return [getattr(self, name) for name in self.names];

    @property
    def elasticity_tensor(self):
        return self.material.elasticity_tensor;

    @property
    def compliance_tensor(self):
        return self.material.compliance_tensor;

