import os.path

import csv
import numpy as np
from HomogenizationWrapper import HomogenizationWrapper

import pyflann

class PatternParameterTable:
    def __init__(self, index_dir, key="compliance", material_file=None):
        self.key = key;
        self.index_dir = index_dir;
        self.material_file = material_file;

        self.lookup_table = self.__load_index(key);
        self.header, self.pattern = self.__load_data("pattern.csv");
        self.material_header, self.materials = self.__load_data("material.csv");

        #self.homogenizer = HomogenizationWrapper(
        #        self.index_dir, self.header, self.material_file);

    def lookup_and_interpolate(self, materials, num_candidates=3,
            rehomogenize=False):
        target_tensors = self.__extract_lookup_keys(materials);
        index, dist = self.lookup_table.nn_index(target_tensors, num_candidates);

        weights = np.ones_like(dist) / dist;
        weights = weights / np.sum(weights, axis=1)[:,np.newaxis];

        param_values = self.pattern[index].astype(float);
        param_values = np.sum(param_values * weights[:,:,np.newaxis], axis=1);

        if rehomogenize:
            assert(self.material_file is not None);
            assert(os.path.exists(self.material_file));
            materials = self.compute_homogenized_materials(param_values);
            material_values = np.array([m.values for m in materials], order="C");
            material_header = ["fitted_{}".format(name)
                    for name in materials[0].names];
        else:
            material_header = ["fitted_{}".format(name)
                    for name in self.material_header];
            material_values = self.materials[index[:,0],:];

        return param_values, material_header, material_values;

    def compute_homogenized_materials(self, param_values):
        raise NotImplementedError("Rehomogenization is out of date");
        materials = [];
        for param in param_values:
            self.homogenizer.homogenize(param);
            materials.append(self.homogenizer.material);
        return materials;

    def __load_index(self, index_name):
        data = self.__load_dataset(index_name);
        table = pyflann.FLANN();
        table.load_index(os.path.join(self.index_dir, index_name + ".index"), data);
        return table;

    def __load_dataset(self, dataset_name):
        data = np.load(os.path.join(self.index_dir, dataset_name + ".npy"));
        return data;

    def __load_data(self, csv_file):
        csv_file = os.path.join(self.index_dir, csv_file);
        data = [];
        with open(csv_file, 'r') as fin:
            reader = csv.reader(fin);
            header = reader.next();
            for row in reader:
                data.append(row);
        return header, np.array(data);

    def __extract_lookup_keys(self, materials):
        if self.key == "compliance":
            target_tensors = np.array([material.compliance_tensor.ravel(order="C")
                for material in materials ]);
        elif self.key == "elasticity":
            target_tensors = np.array([material.elasticity_tensor.ravel(order="C")
                for material in materials ]);
        else:
            raise NotImplementedError("Unknow lookup key: {}".format(self.key));

        return target_tensors;
