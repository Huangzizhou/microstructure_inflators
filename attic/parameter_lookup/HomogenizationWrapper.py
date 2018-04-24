import copy
import datetime
import json
import os
import os.path
import re
from subprocess import check_call, check_output

from microstructures_setting import MICROSTRUCTURES_PATH
from MaterialParameter import MaterialParameter

def update_array(array, param_map):
    for i,val in enumerate(array):
        if isinstance(val, (str, unicode)):
            array[i] = val.format(**param_map);
        elif isinstance(val, list):
            array[i] = update_array(val, param_map);
        elif isinstance(val, dict):
            array[i] = update_dict(val, param_map);
    return array;

def update_dict(config, param_map):
    for key,val in config.iteritems():
        if isinstance(val, (str, unicode)):
            config[key] = val.format(**param_map);
        elif isinstance(val, list):
            config[key] = update_array(val, param_map);
        elif isinstance(val, dict):
            config[key] = update_dict(val, param_map);
    return config;

class HomogenizationWrapper:
    def __init__(self, index_dir, header, material_file):
        self.index_dir = index_dir;
        self.header = header;
        self.material_file = material_file;
        self.__parse_index_path();
        self.__load_modifier_file();

    def __parse_index_path(self):
        full_index_path = os.path.abspath(self.index_dir);
        pattern = "/([23])D/([^/]*)/([^/]*)mm_cell/";
        matcher = re.compile(pattern);
        result = matcher.search(full_index_path);
        assert(result is not None);
        self.dim = int(result.group(1));
        self.wire_file = os.path.join("patterns/{}D/".format(self.dim),
                result.group(2) + ".wire");
        self.cell_size = float(result.group(3));
        #print("dim = {}".format(self.dim));
        #print("wire_file = {}".format(self.wire_file));
        #print("cell_size = {}".format(self.cell_size));

    def __load_modifier_file(self):
        modifier_file = os.path.join(self.index_dir, "lookup.modifier");
        assert(os.path.exists(modifier_file));
        with open(modifier_file) as fin:
            self.modifier_config = json.load(fin);

    def homogenize(self, parameters):
        param_map = {key:val for key,val in zip(self.header, parameters)};
        customized_modifier = copy.deepcopy(self.modifier_config);
        customized_modifier = update_dict(customized_modifier, param_map);

        self.__generate_time_stamp();
        self.__dump_modifier(customized_modifier);
        self.__dump_config();
        self.__tile_and_fit();

    def __generate_time_stamp(self):
        self.tmp_dir = "/tmp";
        self.stamp = datetime.datetime.now().isoformat();

    def __dump_modifier(self, modifier_config):
        self.modifier_file = os.path.join(self.tmp_dir,
                "{}.modifier".format(self.stamp));
        with open(self.modifier_file, 'w') as fout:
            json.dump(modifier_config, fout);

    def __dump_config(self):
        config = {
                "periodic": True,
                "thickness": 0.5,
                "bbox_min": [0.0, 0.0, 0.0][:self.dim],
                "bbox_max": [self.cell_size] * self.dim,
                "wire_network": self.wire_file,
                "repeats": [1, 1, 1][:self.dim],
                "modifier_file": self.modifier_file
                }
        self.config_file = os.path.join(self.tmp_dir,
                "{}.config".format(self.stamp));
        with open(self.config_file, 'w') as fout:
            json.dump(config, fout);

    def __tile_and_fit(self):
        tmp_mesh_file = os.path.join(self.tmp_dir, 
                "{}.msh".format(self.stamp));
        if self.dim == 2:
            exe_name = os.path.join(MICROSTRUCTURES_PATH,
                    "wire_generator/tile_and_fit_periodic_2D.py");
        else:
            exe_name = os.path.join(MICROSTRUCTURES_PATH,
                    "wire_generator/tile_and_fit_periodic.py");
        command = "{} --material {} {} {}".format(exe_name, self.material_file,
                self.config_file, tmp_mesh_file);
        check_output(command.split());

        tmp_json_file = os.path.join(self.tmp_dir,
                "{}_param.json".format(self.stamp));
        assert(os.path.exists(tmp_json_file));
        self.material = MaterialParameter(self.dim, tmp_json_file);

        os.remove(self.config_file);
        os.remove(self.modifier_file);
        if os.path.exists(tmp_mesh_file):
            os.remove(tmp_mesh_file);
        os.remove(tmp_json_file);


