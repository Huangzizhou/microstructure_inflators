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

parser.add_argument('png',    help='png file representing shape')
parser.add_argument('--output', help='file with coefficients for level set representation')
parser.add_argument('--output-png', help='image representing level set representation')
parser.add_argument('--show', action='store_true', help='show image')
parser.add_argument('--output-error', action='store_true', help='output percentage of different pixels')
parser.add_argument('--basisPerDim', type=int, default=5, help='number of points per dimension used')
parser.add_argument('--drop-points', type=float, default=0.0, help='percentage of points to be dropped')


args = parser.parse_args()

img = Image.open(args.png)
arr = np.array(img) / 255 * 2.0 - 1.0

num_rows = len(arr)
num_cols = len(arr[0])

# Make vector flat
data_values = arr.flatten()

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
#print(x1)
#print(x2)
print("Fitting data points with RBF...")
n = args.basisPerDim
if args.drop_points > 0:
    inds = set(random.sample(list(range(len(x1))), int(args.drop_points * len(x1))))

    x1_cheap =          [n for i, n in enumerate(x1)          if i not in inds]
    x2_cheap =          [n for i, n in enumerate(x2)          if i not in inds]
    data_values_cheap = [n for i, n in enumerate(data_values) if i not in inds]

    print("Data points: " + str(len(x1)))
    print("Data points after cut: " + str(len(x1_cheap)))
    rbf = rbf_interpolation.RBFInterpolation(x1_cheap, x2_cheap, data_values_cheap, n, n, n/2)
else:
    rbf = rbf_interpolation.RBFInterpolation(x1, x2, data_values, n, n, n / 2)

# Run function on grid points to obtain a new image array
print("Evaluating points...")
results = rbf(np.array(x1), np.array(x2))
output = []
for i in range(len(x1)):
    if results[i] < 0.0:
        output.append(-1.0)
    else:
        output.append(1.0)

# Transform image array into image
#output_array = np.copy(arr) * 255
output_array = (np.array(output) + 1) / 2.0 * 255
output_flatten = np.array(output)
output_array = output_array.reshape((num_rows, num_cols))
output_array = np.asarray(dtype=np.dtype('uint8'), a=output_array)

# Plot new image
output_img = Image.fromarray(output_array, mode='L')

# Show new image
if args.show:
    output_img.show()

# Save new image
if args.output_png:
    output_img.save(args.output)

# Compute error
if args.output_error:
    error_array = (data_values - output_flatten) / 2
    wrong_pixels = np.sum(np.abs(error_array))
    error = wrong_pixels / len(data_values)
    print("Error percentage: " + str(error))
    print("Wrong pixels: " + str(wrong_pixels))

# Output coefficients
if args.output:
    save(rbf, args.output)




