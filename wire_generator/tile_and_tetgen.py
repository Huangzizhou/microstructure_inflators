#!/usr/bin/env python

import argparse
import hashlib
from subprocess import check_call
import os
import os.path

def generate_stamp(config_file):
    m = hashlib.md5();
    m.update(config_file);
    return m.hexdigest();

def run_tile(config_file, obj_file):
    cmd = "./tile.py -o {} {}".format(obj_file, config_file);
    print(cmd);
    check_call(cmd.split());

def run_tetgen(obj_file, msh_file):
    cmd = "tetgen.py --cmd --flags=\"qpQY\" {} {}".format(obj_file, msh_file);
    print(cmd);
    check_call(cmd.split());

def parse_args():
    parser = argparse.ArgumentParser(
            description="Tile pattern, tetgen it");
    parser.add_argument("config_file", help="configuration file");
    parser.add_argument("msh_file", help="output msh file");
    args = parser.parse_args();
    return args

def main():
    args = parse_args();

    stamp = generate_stamp(args.msh_file);
    tmp_dir = "/tmp"
    tmp_obj = os.path.join(tmp_dir, stamp+".obj");

    run_tile(args.config_file, tmp_obj);
    run_tetgen(tmp_obj, args.msh_file);

    os.remove(tmp_obj);

if __name__ == "__main__":
    main();

