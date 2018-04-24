import LinearElasticitySettings
import PyAssembler

import numpy as np

class MaterialState(object):
    def __init__(self, mesh,
            density=1.0,
            young=5.0,
            poisson=0.0):

        self.mesh = mesh;
        self.density = density;
        self.young_attr_name = "young";
        self.poisson_attr_name = "poisson"

        self.__set_material_attributes(
                np.ones(self.mesh.num_elements) * young,
                np.ones(self.mesh.num_elements) * poisson);
        self.__init_material();

    def __set_material_attributes(self, young, poisson):
        if not self.mesh.has_attribute(self.young_attr_name):
            self.mesh.add_attribute(self.young_attr_name);
        self.mesh.set_attribute(self.young_attr_name, young);

        if not self.mesh.has_attribute(self.poisson_attr_name):
            self.mesh.add_attribute(self.poisson_attr_name);
        self.mesh.set_attribute(self.poisson_attr_name, poisson);

    def __init_material(self):
        self.material = \
                PyAssembler.Material.create_element_wise_isotropic(
                        self.density,
                        self.mesh.raw_mesh,
                        self.young_attr_name,
                        self.poisson_attr_name);

    def update(self, young, poisson):
        young = np.clip(young, 0.1, 5.0);
        poisson = np.clip(poisson, -0.3, 0.3);
        self.__set_material_attributes(young, poisson);
        self.material.update();

    @property
    def material(self):
        return self.__hetero_material;

    @material.setter
    def material(self, mat):
        self.__hetero_material = mat;

    @property
    def young(self):
        return self.mesh.get_attribute(self.young_attr_name).ravel();

    @property
    def poisson(self):
        return self.mesh.get_attribute(self.poisson_attr_name).ravel();

