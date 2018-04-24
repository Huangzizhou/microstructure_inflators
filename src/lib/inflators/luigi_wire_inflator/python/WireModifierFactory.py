import json
import os.path

class WireModifierFactory:
    @classmethod
    def create_from_file(cls, modifier_file):
        config_dir = os.path.dirname(modifier_file);
        with open(modifier_file, 'r') as fin:
            config = json.load(fin);
            config["config_dir"] = config_dir;
            return cls.create_from_dict(config);

    @classmethod
    def create_from_dict(cls, config):
        config_dir = config.get("config_dir", "./");
        modifiers = [];
        if "thickness" in config:
            thickness_config = config["thickness"];
            cls.convert_to_abs_path(thickness_config, config_dir);

            from WireThicknessModifier import WireThicknessModifier
            modifier = WireThicknessModifier(thickness_config);
            modifiers.append(modifier);

        if "vertex_offset" in config:
            offset_config = config["vertex_offset"];
            cls.convert_to_abs_path(offset_config, config_dir);

            from WireVertexOffsetModifier import WireVertexOffsetModifier
            modifier = WireVertexOffsetModifier(offset_config);
            modifiers.append(modifier);
        return modifiers;

    @classmethod
    def convert_to_abs_path(cls, config, config_dir):
        for key,val in config.iteritems():
            if "file" in key and not os.path.isabs(val):
                val = os.path.join(config_dir, val);
                config[key] = val;

