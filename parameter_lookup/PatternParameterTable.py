import os.path

import csv
import numpy as np

import pyflann

class PatternParameterTable:
    def __init__(self, index_dir):
        self.index_dir = index_dir;
        self.elasticity_table = self.__load_index("elasticity");
        self.header, self.pattern = self.__load_data("pattern.csv");

        self.young = self.__load_dataset("young");
        self.poisson = self.__load_dataset("poisson");
        self.shear = self.__load_dataset("shear");

    def lookup(self, materials):
        target_tensors = np.array([material.elasticity_tensor.ravel(order="C")
            for material in materials ]);
        index, dist = self.elasticity_table.nn_index(target_tensors, 1);

        param_values = self.pattern[index];
        young = self.young[index];
        poisson = self.poisson[index];
        shear = self.shear[index];

        return param_values, young, poisson, shear, dist.ravel();

    def lookup_and_interpolate(self, materials):
        target_tensors = np.array([material.elasticity_tensor.ravel(order="C")
            for material in materials ]);
        index, dist = self.elasticity_table.nn_index(target_tensors, 3);

        weights = np.ones_like(dist) / dist;
        weights = weights / np.sum(weights, axis=1)[:,np.newaxis];

        param_values = self.pattern[index].astype(float);
        param_values = np.sum(param_values * weights[:,:,np.newaxis], axis=1);

        index = index[:,0];
        dist = dist[:,0];
        young = self.young[index];
        poisson = self.poisson[index];
        shear = self.shear[index];

        return param_values,  young, poisson, shear, dist.ravel();

    def compute_homogenized_materials(self, param_values):
        # TODO
        pass;

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
