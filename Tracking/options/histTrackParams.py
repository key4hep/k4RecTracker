from itertools import product
from os import getenv
from pathlib import Path
from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import uproot

# general plotting options
labelsize = 20
linewidth = 1.5
majorTickSize = 10
plot_grid_alpha = .7
params = {
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": majorTickSize,  # Tick length
    "ytick.major.size": majorTickSize,
    "xtick.major.width": linewidth,  # Tick line width
    "ytick.major.width": linewidth,
    "axes.linewidth": linewidth,
    "legend.fontsize": labelsize,  # "x-large",
    "axes.labelsize": labelsize,  # "x-large",
    "axes.titlesize": labelsize,  # "x-large",
    "xtick.labelsize": labelsize,  # "x-large",
    "ytick.labelsize": labelsize,  # "x-large",
    "figure.autolayout": True,
}  #'figure.figsize': (15, 5),
mpl.rcParams.update(params)

# argparse
detModOpts = ["v02", "if1", "if2"]
class CaseInsensitiveDict(dict):
    def __getitem__(self, key):
        return super().__getitem__(key.lower())

    def __setitem__(self, key, value):
        super().__setitem__(key.lower(), value)
detModNames = CaseInsensitiveDict({el: el for el in detModOpts})
parser = ArgumentParser()
parser.add_argument(
    "--detectorModel",
    "-m",
    help="Which detector model to run reconstruction for",
    choices=detModOpts + [el.upper() for el in detModOpts],
    type=str,
    default="V02",
)
parser.add_argument(
    "--version", type=str, help="str to identify a run through the pipeline"
)
args = parser.parse_known_args()[0]
corePath = Path(f"{args.version}_{detModNames[args.detectorModel]}")

# strings to build path
processor = "TrackParamExtractor"
basePath = Path(getenv("dtDir", str(Path.home() / "promotion" / "data")))

# Lists to build branch names to be analyzed
trackNames = ["SiTrack", "CluTrack"]
varNames = ["D0", "Omega", "Phi", "TanL", "Z0"]

# build vars based on above vars
keys = [f"{trackName}{varName}" for trackName, varName in product(trackNames, varNames)]
in_file = basePath / processor / corePath.with_suffix(".edm4hep.root")

print(in_file)
with uproot.open(str(in_file) + ":events") as events:
    regex = f"/^({'|'.join(trackNames)})({'|'.join(varNames)})$/"
    data = events.arrays(filter_name=regex, library="pd")
    for var in varNames:
        data[f"d_{var}"] = data[f"{trackNames[0]}{var}"] - data[f"{trackNames[1]}{var}"]

for varName in varNames:
    plt.figure()
    plt.grid(True, which="both", linestyle="--", linewidth=linewidth, alpha=plot_grid_alpha)
    plt.hist(data[f"d_{varName}"], bins=30)
    plt.xlabel(rf"$\Delta$ {varName}")
    plt.ylabel("Frequency")
    plt.title(args.detectorModel)
    plt.show()
