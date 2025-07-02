from argparse import ArgumentParser
from itertools import product
from os import getenv
from pathlib import Path

import awkward as ak
import matplotlib as mpl
import matplotlib.pyplot as plt
import uproot

from commonArgParsing import add_common_args, detModNames, registry

# general plotting options
labelsize = 22
linewidth = 1.5
majorTickSize = 10
plotGridAlpha = 0.7
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

args = add_common_args(ArgumentParser()).parse_known_args()[0]

#################################
# commands to access cov Matrix
#################################
# import ROOT
# ROOT.gInterpreter.LoadFile("edm4hep/utils/cov_matrix_utils.h")
# # ...
# edm4hep.utils.get_cov_value(cov_m, edm4hep.TrackParams.d0, edm4hep.TrackParams.d0)

# Lists to build branch names to be analyzed
trackType = ["SiTrack", "CluTrack"]
varSimilar = ["Phi", "Omega", "TanL"]
varSpread = ["D0", "Z0"]
varNames = varSimilar + varSpread

data = {}

# extract data
for detMod in args.detectorModels:
    corePath = Path(f"{args.version}_{detModNames[detMod]}")

    # strings to build path
    processor = "TrackParamExtractor"
    basePath = Path(getenv("dtDir", str(Path.home() / "promotion" / "data")))

    # build vars based on above vars
    keys = [
        f"{trackName}{varName}" for trackName, varName in product(trackType, varNames)
    ]
    in_file = basePath / processor / corePath.with_suffix(".edm4hep.root")

    with uproot.open(str(in_file) + ":events") as events:
        regex = f"/^({'|'.join(trackType)})({'|'.join(varNames)})$/"
        data[detMod] = events.arrays(filter_name=regex, library="pd")
        for var in varNames:
            data[detMod][f"d_{var}"] = (
                data[detMod][f"{trackType[0]}{var}"]
                - data[detMod][f"{trackType[1]}{var}"]
            )

# # plot data
# for detMod in args.detectorModels:
#     for group, xlim in zip([varSpread, varSimilar], [1.5, 0.0015]):
#         plt.figure()
#         plt.grid(
#             True, which="both", linestyle="--", linewidth=linewidth, alpha=plotGridAlpha
#         )
#         plt.hist(
#             x=[
#                 ak.to_numpy(ak.flatten(data[detMod][f"d_{varName}"]))
#                 for varName in group
#             ],
#             bins=30,
#             label=group,
#             range=(-xlim, xlim),
#         )
#         plt.xlabel(rf"$\Delta$ Si-Clu")
#         plt.ylabel("Frequency")
#         plt.title(f"{detMod}: {','.join(group)}")
#         plt.legend()
#         plt.show()

for type, xlim in zip(trackType, [.03,1]):
    for var in ["D0"]:
        plt.figure()
        plt.grid(
            True, which="both", linestyle="--", linewidth=linewidth, alpha=plotGridAlpha
        )
        plt.hist(
            x=[
                ak.to_numpy(ak.flatten(data[detMod][f"{type}{var}"]))
                for detMod in args.detectorModels
            ],
            bins=30,
            label=[registry.get(detMod).get_name(args) for detMod in args.detectorModels],
            range=(-xlim, xlim),
        )
        plt.ylabel("Frequency")
        plt.title(f"{type}: {var}")
        plt.legend()
        plt.show()
