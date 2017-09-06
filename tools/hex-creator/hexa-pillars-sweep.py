#!/usr/bin/env python
import argparse
import math
import os
import subprocess

import numpy as np

MAX_NUM_PILLARS = 30
MIN_TRIANGLE_SIDE = 0.75
MAX_TRIANGLE_SIDE = 0.99
TRIANGLE_SIDE_EXPERIMENTS = 20
TOLERANCE = 1e-7


def compute_volume_info(experiment_type, p1, p2, p3, p4):
    if experiment_type == "negative":
        height = 1 * math.sqrt(3) / 3.0
        s = 2 * height

        volume_triangle = math.sqrt(3) / 4 * (p1 * s) ** 2
        volume_pillars = p1 * p3 * p4 * s * (1 - math.sqrt(3) / 2 * s * p1)
        volume = volume_triangle + volume_pillars

        total_volume = math.sqrt(3) / 4 * s ** 2

        volume_fraction = volume / total_volume

    else:
        r = p1 * math.sqrt(3) / 3
        x = p1

        volume_triangle = r * p1
        volume_pillars = 2 * math.sqrt(3) / 3 * (1 - p1) * p1 * p3 * p4

        volume_fraction = p1 ** 2 + 2 * (1 - p1) * p1 * p3 * p4

        volume = volume_triangle + volume_pillars
        total_volume = math.sqrt(3) / 3
        assert abs(volume_fraction - volume / total_volume) < TOLERANCE

    print "Volume is: ", volume
    print "Volume fraction is ", volume_fraction

    return volume, volume_fraction


def print_experiment_info(p1, p2, p3, p4):
    print "\nExperimenting with parameters: \n" \
          "p1: " + str(p1) + "\n" \
                             "p2: " + str(p2) + "\n" \
                                                "p3: " + str(p3) + "\n" \
                                                                   "p4: " + str(p4)


def compute_thickness(experiment_type, vol_frac, triangle_side_factor, chirality_factor, num_pillars):
    thickness_ratio = 0.0

    if experiment_type == "negative":
        # compute triangle side
        height = 1 * math.sqrt(3) / 3.0
        s = 2 * height
        triangle_side = triangle_side_factor * s

        # compute volume to be used by pillars
        total_volume = math.sqrt(3) / 3
        volume = total_volume * vol_frac
        volume_pillars = volume - triangle_side ** 2 * math.sqrt(3) / 4  # removes volume already used by hexagon
        if volume_pillars < 0.0:
            print "Warning! Cannot experiment since triangle alone has more volume than the one allowed by our volume fraction."
            return -1.0

        # compute pillar area and thickness
        pillar_area = triangle_side * chirality_factor
        thickness = 2 * volume_pillars / (num_pillars * (2 - triangle_side * math.sqrt(3)))
        total_width = thickness * num_pillars
        thickness_ratio = thickness * num_pillars / pillar_area
        if thickness_ratio > 1.0:
            print "Skipping experiment since pillar area is smaller than necessary width."
            return -1.0

    else:
        r = triangle_side_factor * math.sqrt(3) / 3
        x = triangle_side_factor
        height = math.sqrt(3) / 3 - r

        # compute volume to be used by pillars
        volume_triangle = r * triangle_side_factor
        total_volume = math.sqrt(3) / 3
        volume = total_volume * vol_frac
        volume_pillars = volume - volume_triangle
        if volume_pillars < 0.0:
            print "Warning! Cannot experiment since triangle alone has more volume than the one allowed by our volume fraction."
            return -1.0

        # compute pillar area and thickness
        triangle_side = 2 * x
        pillar_area = triangle_side * chirality_factor

        thickness = volume_pillars / (num_pillars * height)
        thickness_ratio = thickness * num_pillars / pillar_area

        direct_thickness_ratio = (vol_frac - triangle_side_factor ** 2) / (
        2 * (1 - triangle_side_factor) * triangle_side_factor * chirality_factor)
        assert abs(thickness_ratio - direct_thickness_ratio) < TOLERANCE

        if thickness_ratio > 1.0:
            print "Skipping experiment since pillar area is smaller than necessary width."
            return -1.0

    return thickness_ratio


def run_experiment(experiment_type, triangle_side_factor, num_pillars, chirality_factor, thickness_ratio):
    volume, vol_frac = compute_volume_info(experiment_type, triangle_side_factor, num_pillars, chirality_factor, thickness_ratio)

    print "Volume fraction: " + str(vol_frac)

    if experiment_type == "negative":
        run_negative_poisson_experiment(vol_frac, triangle_side_factor, num_pillars, chirality_factor, thickness_ratio)
    else:
        run_positive_poisson_experiment(vol_frac, triangle_side_factor, num_pillars, chirality_factor, thickness_ratio)


def run_negative_poisson_experiment(vol_frac, triangle_side_factor, num_pillars, chirality_factor, thickness_ratio):
    print_experiment_info(triangle_side_factor, num_pillars, chirality_factor, thickness_ratio)
    name = folder_path + '/negative-poisson_volfrac-{}_p1-{}_p2-{}_p3-{}_p4-{}'.format(round(vol_frac, 3),
                                                                            round(triangle_side_factor, 3),
                                                                            round(num_pillars, 3),
                                                                            round(chirality_factor, 3),
                                                                            round(thickness_ratio, 3))

    lock_name = name + '.lock'
    mesh_name = name + '.msh'

    if os.path.isfile(mesh_name):
        print "Already computed"
        pass
    elif os.path.isfile(lock_name):
        print "Locked instance"
    else:
        # lock experiment
        open(lock_name, 'a').close()

        cmd = [cwd + '/auxetic-chiral-creator.py', str(triangle_side_factor), str(num_pillars), str(chirality_factor),
               str(thickness_ratio), name + '.wire', mesh_name]
        subprocess.call(cmd)

        # free experiment
        os.remove(lock_name)


def run_positive_poisson_experiment(vol_frac, triangle_side_factor, num_pillars, chirality_factor, thickness_ratio):
    print_experiment_info(triangle_side_factor, num_pillars, chirality_factor, thickness_ratio)
    name = folder_path + '/positive-poisson_volfrac-{}_p1-{}_p2-{}_p3-{}_p4-{}'.format(round(vol_frac, 3),
                                                                            round(triangle_side_factor, 3),
                                                                            round(num_pillars, 3),
                                                                            round(chirality_factor, 3),
                                                                            round(thickness_ratio, 3))

    lock_name = name + '.lock'
    mesh_name = name + '.msh'

    if os.path.isfile(mesh_name):
        print "Already computed"
        pass
    elif os.path.isfile(lock_name):
        print "Locked instance"
    else:
        # lock experiment
        open(lock_name, 'a').close()

        cmd = [cwd + '/hex-creator-igor-parameters.py', str(triangle_side_factor), str(num_pillars),
               str(chirality_factor),
               str(thickness_ratio), name + '.wire', mesh_name]
        subprocess.call(cmd)

        # free experiment
        os.remove(lock_name)


parser = argparse.ArgumentParser(description='Sweep through different hexa-pillars structures.')
parser.add_argument('--vol-frac', type=float, help='constant volume fraction to be used in experiments')
parser.add_argument('--triangle-side-factor', type=float, help='constant triangle side ratio to be used in experiments')
parser.add_argument('--number-pillars', type=int, help='constant number of pillars to be used in experiments')
parser.add_argument('--chirality-factor', type=float, help='constant chirality factor to be used in experiments')
parser.add_argument('experiment_type', help='negative/positive poisson ratio experiment')
parser.add_argument('instances_folder', help='folder where to store mesh/wire/param files')
parser.add_argument('table_file', help='table with elasticity properties from experimented microstructures')

args = parser.parse_args()
cwd = os.getcwd()
folder_path = args.instances_folder

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

experiment_type = args.experiment_type


def test(experiment_type, vol_frac, triangle_side_factor, number_pillars, chirality_factor):
    thickness_ratio = compute_thickness(experiment_type, vol_frac, triangle_side_factor, chirality_factor,
                                        number_pillars)
    if thickness_ratio >= 0.0:
        run_experiment(experiment_type, triangle_side_factor, number_pillars, chirality_factor, thickness_ratio)

        compute_volume_info(experiment_type, triangle_side_factor, number_pillars, chirality_factor, thickness_ratio)


# First, verify if volume fraction is set
if not args.vol_frac is None:
    vol_frac = args.vol_frac

    if not args.chirality_factor is None:
        chirality_factor = args.chirality_factor

        if not args.triangle_side_factor is None:
            triangle_side_factor = args.triangle_side_factor

            # experiment with different numbers of pillars
            num_pillar_values = range(1, MAX_NUM_PILLARS + 1, 1)
            for index, number_pillars in enumerate(num_pillar_values):
                test(experiment_type, vol_frac, triangle_side_factor, number_pillars, chirality_factor)

        else:

            if not args.number_pillars is None:
                number_pillars = args.number_pillars

                # experiment with different sizes of triangles
                triangle_side_values = np.linspace(MIN_TRIANGLE_SIDE, MAX_TRIANGLE_SIDE, TRIANGLE_SIDE_EXPERIMENTS)
                for index, triangle_side_factor in enumerate(triangle_side_values):
                    test(experiment_type, vol_frac, triangle_side_factor, number_pillars, chirality_factor)

            else:
                print "Warning: testing in 2 dimensions with triangle sides and pillars"

                triangle_side_values = np.linspace(MIN_TRIANGLE_SIDE, MAX_TRIANGLE_SIDE, TRIANGLE_SIDE_EXPERIMENTS)
                for index, triangle_side_factor in enumerate(triangle_side_values):
                    num_pillar_values = range(1, MAX_NUM_PILLARS + 1, 1)
                    for index, number_pillars in enumerate(num_pillar_values):
                        test(experiment_type, vol_frac, triangle_side_factor, number_pillars, chirality_factor)


    else:
        print('Warning! Currently, it is necessary to provide chirality factor for experiments.')

else:
    if not args.chirality_factor is None:
        chirality_factor = args.chirality_factor
        print "Warning: testing in 3 dimensions with volfrac, triangle sides and pillars"

        num_pillar_values = range(10, 50, 10)
        for index, number_pillars in enumerate(num_pillar_values):
            triangle_side_values = [0.8, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98]
            for index, triangle_side_factor in enumerate(triangle_side_values):
                vol_frac_values = [0.9, 0.92, 0.94, 0.96, 0.98]
                for vol_frac in vol_frac_values:
                    test(experiment_type, vol_frac, triangle_side_factor, number_pillars, chirality_factor)
    else:
        chirality_factor_values = [0.6, 0.7, 0.8, 0.9]
        for chirality_factor in chirality_factor_values:
            print "Warning: testing in 4 dimensions with volfrac, triangle sides, chirality and pillars"

            num_pillar_values = range(10, 50, 10)
            for index, number_pillars in enumerate(num_pillar_values):
                triangle_side_values = [0.82, 0.9, 0.98]
                for index, triangle_side_factor in enumerate(triangle_side_values):
                    vol_frac_values = [0.9, 0.92, 0.94, 0.96, 0.98]
                    for vol_frac in vol_frac_values:
                        test(experiment_type, vol_frac, triangle_side_factor, number_pillars, chirality_factor)
