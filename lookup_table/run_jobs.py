#!/usr/bin/env python

import argparse
import json
import os
import os.path
from subprocess import check_call

def load_job_setting(job_file):
    with open(job_file, 'r') as fin:
        job_setting = json.load(fin);
        return job_setting;

def clean(out_dir):
    command = "rm -rf {}".format(out_dir);
    print(command);
    check_call(command.split());
    command = "mkdir -p {}".format(out_dir);
    print(command);
    check_call(command.split());

def run_job(job_file):
    command = "submit_job.py --batch-size=3 {}".format(job_file);
    print(command);
    check_call(command.split());

def generate_index(config_dir, out_dir, index_dir):
    microstructure_path = os.environ["MICROSTRUCTURES_PATH"];
    command = "{}/parameter_lookup/generate_lookup_table.py --dim 2 -o {} {} {}".format(
            microstructure_path, index_dir, config_dir, out_dir);
    print(command);
    check_call(command.split());

def get_out_dir(setting, job_dir):
    out_dir = os.path.join(setting["root_dir"], setting["out_dir"]);
    if not os.path.isabs(out_dir):
        out_dir = os.path.join(job_dir, out_dir);
    return out_dir;

def parse_args():
    parser = argparse.ArgumentParser(
            description="Process and run given job script");
    parser.add_argument("action", type=str, choices=("run", "clean", "index"));
    parser.add_argument("job_files", nargs="+");
    args = parser.parse_args();
    return args;

def main():
    args = parse_args();
    for job_file in args.job_files:
        assert(os.path.exists(job_file));
        setting = load_job_setting(job_file);
        out_dir = get_out_dir(setting, os.path.dirname(job_file));
        base_dir = os.path.abspath(os.path.join(out_dir, "../"));

        if args.action == "clean":
            clean(out_dir);
        elif args.action == "run":
            run_job(job_file);
        elif args.action == "index":
            config_dir = os.path.join(base_dir, "configs");
            index_dir = os.path.join(base_dir, "index");
            clean(index_dir);
            generate_index(config_dir, out_dir, index_dir);
        else:
            raise NotImplementedError("Unknown action: {}".format(args.action));

if __name__ == "__main__":
    main();
