#!/usr/bin/env python3
import argparse
import random

import numpy as np
from PIL import Image

import rbf_interpolation


def save(rbf, output_path):
    out = open(output_path, 'w')

    # Method
    out.write("rbf\n")

    # Properties of basis functions
    out.write("{}\n".format(rbf.epsilon))

    # Coefficients
    coeffs = rbf.coeffs
    out.write("{}\t{}\n".format(rbf.d1, rbf.d2))

    coeffs_list = list(coeffs)
    coeffs_string = ""
    for i in range(0, len(coeffs_list)):
        coeffs_string += '{:.16f}'.format(coeffs_list[i]) + "\t"

    out.write(coeffs_string)

    out.close()


parser = argparse.ArgumentParser(description='Transform png into level set representation')

parser.add_argument('frame1',    help='png file representing initial shape')
parser.add_argument('frame2',    help='png file representing final shape')
parser.add_argument('--output-prefix', default="tmp_", help='prefix for output')
parser.add_argument('--frames-number', type=int, default=5, help='number of intermediate frames')
parser.add_argument('--basis-per-dim', type=int, default=3, help='number of points per dimension used')

args = parser.parse_args()

img1 = Image.open(args.frame1)
img2 = Image.open(args.frame2)
arr1 = np.array(img1) / 255 * 2.0 - 1.0
arr2 = np.array(img2) / 255 * 2.0 - 1.0

num_rows = len(arr1)
num_cols = len(arr1[0])

# Make vector flat
frame1_values = arr1.flatten()
frame2_values = arr2.flatten()

# Construct grid with positions between -1 to 1 on both directions
# Given pixel size, first data point should be at -1 + pixel_size/2. Last should be at 1 - pixel_size/2. All the other N - 2
# in between
pixel_size1 = (1.0 - (-1.0)) / num_cols
pixel_size2 = (1.0 - (-1.0)) / num_rows
x1 = []
x2 = []
for i in range(num_rows):
    x2_pos = (1.0 - pixel_size2 / 2.0) - i * pixel_size2
    for j in range(num_cols):
        x1_pos = (-1.0 + pixel_size1/2.0) + j * pixel_size1

        x1.append(x1_pos)
        x2.append(x2_pos)

# Use RBF interpolation class to build a function approximating our shape
n = args.basis_per_dim
print("Fitting data points of first frame with RBF...")
rbf1 = rbf_interpolation.RBFInterpolation(x1, x2, frame1_values, n, n, n / 2)
coeffs1 = np.array(rbf1.coeffs)
print("Fitting data points of second frame with RBF...")
rbf2 = rbf_interpolation.RBFInterpolation(x1, x2, frame2_values, n, n, n / 2)
coeffs2 = np.array(rbf2.coeffs)

# Now, for each intermediate frame, produces linear interpolation
t_values = np.linspace(0, 1, args.frames_number)
for t in t_values:

    # for each coefficient, run linear interpolation
    intermediate_coeffs = t * coeffs1 + (1.0 - t) * coeffs2
    print(intermediate_coeffs)

    rbf = rbf_interpolation.RBFInterpolation(d1=n, d2=n, epsilon=n / 2, coeffs=list(intermediate_coeffs))

    results = rbf(np.array(x1), np.array(x2))
    output = []
    print(results)
    for i in range(len(x1)):
        if results[i] < 0.0:
            output.append(-1.0)
        else:
            output.append(1.0)

    # Transform image array into image
    output_array = (np.array(output) + 1) / 2.0 * 255
    output_flatten = np.array(output)
    output_array = output_array.reshape((num_rows, num_cols))
    output_array = np.asarray(dtype=np.dtype('uint8'), a=output_array)

    output_img = Image.fromarray(output_array, mode='L')
    output_img.save(args.output_prefix + str(t) + ".png")
