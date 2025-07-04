from argparse import ArgumentParser
from itertools import product
from os import getenv
from pathlib import Path

import awkward as ak
import matplotlib as mpl
import matplotlib.pyplot as plt
import uproot

from commonArgParsing import add_common_args, detModNames, registry
from plotting import linewidth, params, plotGridAlpha

plt.style.use("seaborn-v0_8-colorblind")

mpl.rcParams.update(params)

parser = add_common_args(ArgumentParser())
plot_opts = parser.add_argument_group("Plotting opts", "which plots should be shown")
plot_opts.add_argument(
    "--track",
    action="store_true",
    help="Show difference between Silicon and Clupatra tracks",
)
plot_opts.add_argument(
    "--detmods", action="store_true", help="Show difference between detector models"
)
args = parser.parse_known_args()[0]

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

if args.track:
    # plot data
    for detMod in args.detectorModels:
        for group, xlim in zip([varSpread, varSimilar], [1.5, 0.0015]):
            plt.figure()
            plt.grid(
                True,
                which="both",
                linestyle="--",
                linewidth=linewidth,
                alpha=plotGridAlpha,
            )
            plt.hist(
                x=[
                    ak.to_numpy(ak.flatten(data[detMod][f"d_{varName}"]))
                    for varName in group
                ],
                bins=30,
                label=group,
                range=(-xlim, xlim),
            )
            plt.xlabel(rf"$\Delta$ Si-Clu")
            plt.ylabel("Frequency")
            plt.title(
                rf"Diff Si-Clu in $\mathtt{{{registry.get(detMod).get_name(args)}}}$: {','.join(group)}"
            )
            plt.legend()
            plt.show()

if args.detmods:
    xlims = {type: None for type in trackType}
    xlims["SiTrack"] = {"D0": 0.03, "Omega": 0.0002}
    xlims["CluTrack"] = {"D0": 0.8, "Omega": 0.00025}
    for type in trackType:
        for var in ["D0", "Omega"]:
            plt.figure()
            plt.grid(
                True,
                which="both",
                linestyle="--",
                linewidth=linewidth,
                alpha=plotGridAlpha,
            )
            plt.hist(
                x=[
                    ak.to_numpy(ak.flatten(data[detMod][f"{type}{var}"]))
                    for detMod in args.detectorModels
                ],
                bins=30,
                label=[
                    rf"$\mathtt{{{registry.get(detMod).get_name(args)}}}$"
                    for detMod in args.detectorModels
                ],
                range=(-xlims[type][var], xlims[type][var])
                if var in xlims[type]
                else None,
            )
            plt.ylabel("Frequency")
            plt.title(f"Diff DetMods $\mathtt{{{type}}}$: {var}")
            plt.legend()
            plt.show()
