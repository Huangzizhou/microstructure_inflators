import math
import numpy as np
from WireAttribute import WireAttribute
from utils.quaternion import Quaternion

import PyMesh

class WireSymmetryAttribute(WireAttribute):
    def __init__(self):
        super(WireSymmetryAttribute, self).__init__();
        self.symmetry_type = "orthotropic";

    def set_symmetry_type(self, symmetry_type):
        if (symmetry_type == "orthotropic"):
            self.symmetry_type = symmetry_type;
        elif (symmetry_type == "isotropic"):
            self.symmetry_type = symmetry_type;
        else:
            raise NotImplementedError(
                    "Symmetry type ({}) is not supported".format(
                        symmetry_type));

    def compute(self, wire_network):
        self.wire_network = wire_network;

        bbox_min, bbox_max = wire_network.bbox;
        bbox_center = 0.5 * (bbox_min + bbox_max);
        bbox_size = bbox_max - bbox_min;

        self.__initialize_grid(0.001 * np.amax(bbox_size));
        if self.wire_network.dim == 2:
            if (self.symmetry_type == "orthotropic"):
                self.__initialize_reflective_symmetries_2D(bbox_center);
            elif (self.symmetry_type == "isotropic"):
                self.__initialize_isotropic_symmetries_2D(bbox_center);
            else:
                assert(False);
        elif self.wire_network.dim == 3:
            if (self.symmetry_type == "orthotropic"):
                self.__initialize_reflective_symmetries_3D(bbox_center);
            elif (self.symmetry_type == "isotropic"):
                self.__initialize_isotropic_symmetries_3D(bbox_center);
            else:
                assert(False);
        else:
            raise NotImplementedError("Unknown dimension: {}".format(
                self.wire_network.dim));

    def __initialize_grid(self, cell_size):
        self.grid = PyMesh.HashGrid.create(cell_size, self.wire_network.dim);
        for i,v in enumerate(self.wire_network.vertices):
            self.grid.insert(i, v);

    def __initialize_reflective_symmetries_2D(self, bbox_center):
        X = np.array([-1, 1]);
        Y = np.array([ 1,-1]);
        #XY= np.array([-1,-1]);

        self.symmetries = [
                lambda v: (v-bbox_center)*X + bbox_center,
                lambda v: (v-bbox_center)*Y + bbox_center,
                #lambda v: (v-bbox_center)*XY + bbox_center
                ];

    def __initialize_isotropic_symmetries_2D(self, bbox_center):
        X = np.array([-1, 1]);
        Y = np.array([ 1,-1]);

        rot_mat_gen = lambda (angle) : np.array([
            [math.cos(angle), -math.sin(angle)],
            [math.sin(angle),  math.cos(angle)] ]);

        rot_90 = rot_mat_gen(math.pi * 0.5);
        rot_180 = rot_mat_gen(math.pi);

        self.symmetries = [
                lambda v: (v-bbox_center)*X + bbox_center,
                lambda v: (v-bbox_center)*Y + bbox_center,
                lambda v: rot_90.dot(v-bbox_center) + bbox_center,
                lambda v: rot_180.dot(v-bbox_center) + bbox_center,
                ];

    def __initialize_reflective_symmetries_3D(self, bbox_center):
        X = np.array([-1, 1, 1]);
        Y = np.array([ 1,-1, 1]);
        Z = np.array([ 1, 1,-1]);
        XY= np.array([-1,-1, 1]);
        YZ= np.array([ 1,-1,-1]);
        ZX= np.array([-1, 1,-1]);
        XYZ=np.array([-1,-1,-1]);

        self.symmetries = [
                lambda v: (v-bbox_center)*X + bbox_center,
                lambda v: (v-bbox_center)*Y + bbox_center,
                lambda v: (v-bbox_center)*Z + bbox_center,
                lambda v: (v-bbox_center)*XY + bbox_center,
                lambda v: (v-bbox_center)*YZ + bbox_center,
                lambda v: (v-bbox_center)*ZX + bbox_center,
                lambda v: (v-bbox_center)*XYZ + bbox_center,
                ];

    def __initialize_isotropic_symmetries_3D(self, bbox_center):
        """ Captures all symmetries processed by a cube.
        """
        X = np.array([-1, 1, 1]);
        Y = np.array([ 1,-1, 1]);
        Z = np.array([ 1, 1,-1]);

        rot_X_90 = Quaternion.fromAxisAngle(
                [1.0, 0.0, 0.0], math.pi * 0.5).to_matrix();
        rot_Y_90 = Quaternion.fromAxisAngle(
                [0.0, 1.0, 0.0], math.pi * 0.5).to_matrix();
        rot_Z_90 = Quaternion.fromAxisAngle(
                [0.0, 0.0, 1.0], math.pi * 0.5).to_matrix();

        rot_X_180 = Quaternion.fromAxisAngle(
                [1.0, 0.0, 0.0], math.pi).to_matrix();
        rot_Y_180 = Quaternion.fromAxisAngle(
                [0.0, 1.0, 0.0], math.pi).to_matrix();
        rot_Z_180 = Quaternion.fromAxisAngle(
                [0.0, 0.0, 1.0], math.pi).to_matrix();

        self.symmetries = [
                lambda v: (v-bbox_center)*X + bbox_center,
                lambda v: (v-bbox_center)*Y + bbox_center,
                lambda v: (v-bbox_center)*Z + bbox_center,
                lambda v: rot_X_90.dot(v-bbox_center) + bbox_center,
                lambda v: rot_Y_90.dot(v-bbox_center) + bbox_center,
                lambda v: rot_Z_90.dot(v-bbox_center) + bbox_center,
                lambda v: rot_X_180.dot(v-bbox_center) + bbox_center,
                lambda v: rot_Y_180.dot(v-bbox_center) + bbox_center,
                lambda v: rot_Z_180.dot(v-bbox_center) + bbox_center,
                ];

    @property
    def value(self):
        raise NotImplementedError("This attribute is abstract");
