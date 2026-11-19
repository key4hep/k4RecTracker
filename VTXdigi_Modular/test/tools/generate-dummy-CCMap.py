import argparse

import numpy as np

NORMALISATION = 1.0
# N_VOXELS and MATRIX_SIZE cannot be changed, the logic of GetMatrix depends on these exact numbers
N_VOXELS = [3, 3, 3]  # u v w
MATRIX_SIZE = 5


def GetMatrix(i_u, i_v, i_w, matrix_size, normalisation, printMatrix=False):
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


def WriteHeader(file, pixel_dimensions, n_voxels, normalisation):
    file.write("this is just a header that will be filled\n")
    file.write(f"internal {normalisation}\n")
    file.write("##TURN## ##TILT## 1.0\n")
    file.write("0.0 0.0 0.0\n")
    file.write(
        f"{pixel_dimensions[0]} {pixel_dimensions[1]} {pixel_dimensions[2]} 0.0 0.0 0.0 0.0 {n_voxels[0]:n} {n_voxels[1]:n} {n_voxels[2]:n} 0.0\n"
    )
    return


def main():
    parser = argparse.ArgumentParser(
        description="Generate a dummy charge collection map for VTXdigi_Modular tests",
    )
    parser.add_argument("--fileName", help="output file", required=True)
    parser.add_argument(
        "--pitch", type=float, metavar="UM", help="pixel pitch in um", required=True
    )
    parser.add_argument(
        "--thickness",
        type=float,
        metavar="UM",
        help="sensor thickness in um",
        required=True,
    )
    args = parser.parse_args()
    print(f"[generate-dummy-CCmap.py] Generating file {args.fileName}")
    print(
        f"[generate-dummy-CCmap.py] with pitch {args.pitch} um and thickness {args.thickness} um"
    )

    pixel_dimensions = [args.thickness, args.pitch, args.pitch]  # thickness, pitch_u, pitch_v

    with open(args.fileName, "w") as file:
        WriteHeader(file, pixel_dimensions, N_VOXELS, NORMALISATION)

        for i_u in range(N_VOXELS[0]):
            for i_v in range(N_VOXELS[1]):
                for i_w in range(N_VOXELS[2]):
                    # line begins with indices, CCMap is 1-indexed...
                    line = f"{i_u + 1:n} {i_v + 1:n} {i_w + 1:n}"

                    matrix = GetMatrix(i_u, i_v, i_w, MATRIX_SIZE, NORMALISATION)

                    # matrix is stored row-wise in the line
                    for j_v in range(MATRIX_SIZE):
                        line += " " + " ".join(f"{val:.4f}" for val in matrix[:, j_v])
                    line += "\n"
                    _ = file.write(line)


if __name__ == "__main__":
    main()
