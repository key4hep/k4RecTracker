import argparse
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")  # block windows from opening
import matplotlib.pyplot as plt

HEADER_LINES = 5


def ReadCCMap(fileName):
    """Read a CCMap in Allpix Squared format. Returns the sensor thickness (um), the voxel counts
    (u, v, w) and an array of shape (n_u, n_v, n_w) holding the sum of each voxel's matrix."""
    with open(fileName) as file:
        lines = file.readlines()

    header = lines[HEADER_LINES - 1].split()
    if len(header) != 11:
        raise ValueError(f"5th header line has {len(header)} entries, expected 11")
    thickness = float(header[0])
    n_voxels = [int(header[7]), int(header[8]), int(header[9])]

    matrixSums = np.full(n_voxels, np.nan)
    for i_line, line in enumerate(lines[HEADER_LINES:], start=HEADER_LINES + 1):
        if not line.strip() or line.startswith("#"):
            continue
        entries = line.split()
        # CCMap voxel indices are 1-indexed
        i_u, i_v, i_w = (int(entry) - 1 for entry in entries[:3])
        if not np.isnan(matrixSums[i_u, i_v, i_w]):
            raise ValueError(f"duplicate voxel ({i_u}, {i_v}, {i_w}) in line {i_line}")
        matrixSums[i_u, i_v, i_w] = np.sum(np.array(entries[3:], dtype=float))

    if np.isnan(matrixSums).any():
        missing = np.argwhere(np.isnan(matrixSums))
        raise ValueError(f"{len(missing)} voxels missing in CCMap, eg. {tuple(missing[0])}")

    return thickness, n_voxels, matrixSums


def PlotProfile(outputFile, title, w_edges, profile, meanCollectionDepth):
    fig, ax = plt.subplots(figsize=(7, 4.5))

    ax.stairs(profile, w_edges, color="C0", label="CCMap")
    ax.stairs(profile, w_edges, color="C0", alpha=0.15, fill=True)

    ax.axvline(
        0.0,
        color="gray",
        linestyle=":",
    )
    ax.axvline(meanCollectionDepth, color="C1", label="Mean collection depth")

    ax.set_xlim(w_edges[0], w_edges[-1])
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("w (sensor local) [um]")
    ax.set_ylabel("collected charge / deposited charge")
    ax.set_title(title)
    ax.grid()
    ax.legend()

    fig.savefig(outputFile)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", help="input CCMap file", required=True)
    parser.add_argument(
        "--output", "-o", help="output plot (any matplotlib format).", required=True
    )
    args = parser.parse_args()

    thickness, n_voxels, matrixSums = ReadCCMap(args.input)

    # w = 0 at the sensor centre, as in VTXdigi_Modular
    w_edges = np.linspace(-0.5 * thickness, 0.5 * thickness, n_voxels[2] + 1)
    w_centers = 0.5 * (w_edges[:-1] + w_edges[1:])

    # per w-slice: average of the matrix sums over all (u, v) voxels in that slice
    profile = matrixSums.mean(axis=(0, 1))
    if profile.sum() <= 0.0:
        raise ValueError("CCMap collects no charge (all matrices sum to 0)")

    meanCollectionDepth = np.sum(profile * w_centers) / np.sum(profile)
    keptAfterShift = np.abs(w_centers - meanCollectionDepth) < 0.5 * thickness
    chargeLoss = 1.0 - profile[keptAfterShift].sum() / profile.sum()

    print(f"  CCMap                       {args.input}")
    print(f"  thickness                   {thickness} um")
    print(f"  voxel count (u, v, w)       {n_voxels[0]}, {n_voxels[1]}, {n_voxels[2]}")
    print(f"  mean collection depth       {meanCollectionDepth:.3g} um")
    print(f"  charge loss when shifting   {chargeLoss * 100:.3g} %")

    PlotProfile(args.output, Path(args.input).name, w_edges, profile, meanCollectionDepth)
    print(f"Saved plot to                 {args.output}")


if __name__ == "__main__":
    main()
