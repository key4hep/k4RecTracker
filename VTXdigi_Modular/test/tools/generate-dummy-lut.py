import os

import numpy as np
import pandas as pd
from numpy import array as arr

output_file = "lookup_table_dummy.init"

n_voxels = [3, 3, 3]  # u v w
matrix_size = 5
pixel_size = [43.5, 25.0, 25.0]  # thickness, pitch_u, pitch_v
normalisation = 1.0


def GetMatrix(i_u, i_v, i_w, printMatrix=False):
    w_efficiencies = [
        0.01,
        0.2,
        0.8,
    ]  # sensor gets inefficient the further the charge is deposited from the surface
    efficiency = w_efficiencies[i_w]

    matrix = np.zeros([matrix_size, matrix_size])

    # for the 3x3 voxels in u/v we define 3 cases
    # corner | edge   | corner
    # edge   | center | edge
    # corner | edge   | corner
    if i_u == 1 and i_v == 1:  # center
        matrix[2, 2] = 1.0
    elif i_u == 1 or i_v == 1:  # edge
        matrix[2, 2] = 1.0

        ## find out on which edge we are
        if i_u == 1:
            if i_v == 0:
                matrix[2, 1] = 0.5  # left
            else:
                matrix[2, 3] = 0.5  # right
        else:
            if i_u == 0:
                matrix[1, 2] = 0.5  # down
            else:
                matrix[3, 2] = 0.5  # up
    else:  # corner
        matrix[2, 2] = 1.0

        # direct neighbor sharing
        if i_u == 0:
            matrix[1, 2] = 0.5
        elif i_u == 2:
            matrix[3, 2] = 0.5

        if i_v == 0:
            matrix[2, 1] = 0.5
        elif i_v == 2:
            matrix[2, 3] = 0.5

        # sharing across corner
        if i_u == 0 and i_v == 0:
            matrix[1, 1] = 0.25
        elif i_u == 0 and i_v == 2:
            matrix[1, 3] = 0.25
        elif i_u == 2 and i_v == 0:
            matrix[3, 1] = 0.25
        elif i_u == 2 and i_v == 2:
            matrix[3, 3] = 0.25

    matrix_sum = np.sum(matrix)
    matrix = matrix / matrix_sum * efficiency * normalisation
    if printMatrix:
        print("j_v | matrix")
        for j_v in range(matrix_size):
            j_v = matrix_size - 1 - j_v
            print(" ", j_v, "|" + " ".join(f"{val:6.3f}" for val in matrix[:, j_v]))
    return matrix


def WriteHeader(file):
    f.write(f"this is just a header that will be filled\n")
    f.write(f"internal {normalisation}\n")
    f.write(f"##TURN## ##TILT## 1.0\n")
    f.write(f"0.0 0.0 0.0\n")
    f.write(
        f"{pixel_size[0]} {pixel_size[1]} {pixel_size[2]} 0.0 0.0 0.0 0.0 {n_voxels[0]:n} {n_voxels[1]:n} {n_voxels[2]:n} 0.0\n"
    )
    return


with open(output_file, "w") as f:
    WriteHeader(f)

    for i_u in range(n_voxels[0]):
        for i_v in range(n_voxels[1]):
            for i_w in range(n_voxels[2]):
                # line begins with indices, LUT is 1-indexed...
                line = f"{i_u + 1:n} {i_v + 1:n} {i_w + 1:n}"

                matrix = GetMatrix(i_u, i_v, i_w)
                # matrix is stored row-wise in the line
                for j_v in range(matrix_size):
                    line += " " + " ".join(f"{val:.4f}" for val in matrix[:, j_v])
                line += "\n"
                _ = f.write(line)
