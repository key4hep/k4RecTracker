#############################################
# call with `python3` NOT `k4run`
#############################################

from argparse import ArgumentParser
from itertools import product
from os import getenv
from pathlib import Path

import awkward as ak
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import uproot

from commonArgParsing import add_common_args, detModNames, registry
from plotting import my_line_styles
from utils import is_outlier

plt.style.use(["seaborn-v0_8-colorblind", "vics_basic"])

threshold_outlier_detection = 4

my_hist_type = "step"
my_line_width = 2.5
my_n_bins = 30

#############################################
# arg parsing
#############################################

# import common args
parser = add_common_args(ArgumentParser())
parser.add_argument(
    "--mode",
    choices=["nominal", "uncertainty", "both"],
    default="nominal",
    help="Choose what to process: 'nominal' (default), 'uncertainty', or 'both'.",
)
parser.add_argument("--debug", action="store_true")
parser.add_argument("--rm-outliers", action="store_true")
# add plotting options
plot_opts = parser.add_argument_group("Plotting opts", "which plots should be shown")
plot_opts.add_argument(
    "--track",
    action="store_true",
    help="Show difference between Silicon and Clupatra tracks",
)
plot_opts.add_argument(
    "--detmods", action="store_true", help="Show difference between detector models"
)

# parse args
args = parser.parse_known_args()[0]
assert len(parser.parse_known_args()[1]) == 0, (
    f"Unknown args provided: {parser.parse_known_args()[1]}"
)


#############################################
# Lists to build branch names to be analyzed
#############################################

track_types = ["SiTrack", "CluTrack"]
var_similar = ["Phi", "Omega", "TanL"]
var_spread = ["D0", "Z0"]
var_names = var_similar + var_spread

#############################################
# extract data from root file
#############################################
data = {}
for detMod in args.detectorModels:
    corePath = Path(f"{args.version}_{detModNames[detMod]}")

    # strings to build path
    processor = "TrackParamExtractor"
    basePath = Path(getenv("dtDir", str(Path.home() / "promotion" / "data")))

    # build vars based on above vars
    keys = [f"{trackName}{varName}" for trackName, varName in product(track_types, var_names)]
    in_file = basePath / processor / corePath.with_suffix(".edm4hep.root")

    with uproot.open(str(in_file) + ":events") as events:
        # regex to match desired branch names
        regex = (
            f"/^({'|'.join(track_types)}){'(Unc)?' if args.mode != 'nominal' else ''}"
            + f"({'|'.join(var_names)})$/"
        )
        if args.debug:
            print(f"Regex to match branches we are interested in: {regex}")
        data[detMod] = events.arrays(filter_name=regex, library="pd")
        for var in var_names:
            data[detMod][f"d_{var}"] = (
                data[detMod][f"{track_types[0]}{var}"] - data[detMod][f"{track_types[1]}{var}"]
            )
        if args.debug:
            print("Matched branches are:")
            print(data[detMod].columns)


#############################################
# plotting funcs
#############################################
def process_data_for_hist(data, det_mod, var_column_name, rm_outliers, thresh_outlier_detection):
    # TODO: remove comment?
    """
    Process the data for plotting, handling outliers if required.

    Args:
        data: The full data dictionary.
        det_mod: The detector model name.
        var_column_name: The variable name.
        rm_outliers: Boolean flag to remove outliers.
        thresh_outlier_detection: Threshold for outlier detection.

    Returns:
        numpy array: Processed data ready for plotting.
    """
    data_array = ak.to_numpy(ak.flatten(data[det_mod][f"{var_column_name}"]))

    if rm_outliers:
        # Remove outliers
        return data_array[~is_outlier(data_array, thresh=thresh_outlier_detection)]

    return data_array


def debug_info_outlier_removal(
    data, thresh, args, var_list, fixed_data_key, det_mod_variable: bool
):
    if args.rm_outliers and args.debug:
        star_string = "*" * 50
        print(f"\n\n{star_string}\nOutlier threshold value is: {thresh}\n{star_string}")
        for current_var in var_list:
            data_to_clean = (
                data[current_var][fixed_data_key]
                if det_mod_variable
                else data[fixed_data_key][current_var]
            )
            pre_processed_data = ak.to_numpy(ak.flatten(data_to_clean))
            print(
                f"{np.sum(is_outlier(pre_processed_data))} outliers will be removed out of {len(pre_processed_data)} values for {current_var}"
            )


# general plotting func
def plot_track_param_hist(data, thresh_outlier_detection, args, var, labels, hist_args):
    if var:  # TODO: support for var_group case
        debug_info_outlier_removal(
            data,
            thresh_outlier_detection,
            args,
            args.detectorModels,
            var,
            det_mod_variable=True,
        )
    with mpl.rc_context({"axes.titlesize": 18}):
        plt.figure()
        plt.hist(
            **hist_args,
            bins=my_n_bins,
            histtype=my_hist_type,
            linewidth=my_line_width,
            linestyle=my_line_styles,
        )
        plt.xlabel(labels["xlabel"] if "xlabel" in labels else None)
        plt.ylabel("Frequency")
        plt.suptitle(labels["suptitle"])
        plt.title(labels["title"])
        plt.legend()
        plt.show()


# plotting func for collective plot of group of vars
def plot_track_param_hist_var_groups(data, thresh_outlier_detection, args, det_mod, labels, group):
    hist_args_diff_detmods = {
        "x": [
            process_data_for_hist(
                data,
                det_mod,
                f"d_{var_name}",
                args.rm_outliers,
                thresh_outlier_detection,
            )
            for var_name in group
        ],
        "label": group,
    }
    plot_track_param_hist(
        data,
        thresh_outlier_detection,
        args,
        None,
        labels,
        hist_args_diff_detmods,
    )


# plotting func for collective plot of all det mods
def plot_track_param_hist_diff_detmods(data, thresh_outlier_detection, args, var, labels):
    hist_args_diff_detmods = {
        "x": [
            process_data_for_hist(data, det_mod, var, args.rm_outliers, thresh_outlier_detection)
            for det_mod in args.detectorModels
        ],
        "label": [registry.get(det_mod).get_name(args.detname) for det_mod in args.detectorModels],
    }
    plot_track_param_hist(
        data,
        thresh_outlier_detection,
        args,
        var,
        labels,
        hist_args_diff_detmods,
    )


#############################################
# actual plotting (calling the funcs)
#############################################

if args.mode != "nominal":  # uncertainty or both
    for track_type in [track_types[0]]:
        for var in var_names:
            # define labels
            uncertainty_labels = {
                "xlabel": r"$\sigma$",
                "suptitle": rf"$\sigma$({var})",
                "title": rf"single $\mu$, {track_type} {' (no outliers)' if args.rm_outliers else ''}",
            }
            # plotting ;)
            plot_track_param_hist_diff_detmods(
                data,
                threshold_outlier_detection,
                args,
                f"{track_type}Unc{var}",
                uncertainty_labels,
            )


if args.mode != "uncertainty":  # nominal or both
    if args.track:
        for det_mod in args.detectorModels:
            for group in [var_spread, var_similar]:
                # define labels
                nominal_var_groups_labels = {
                    "xlabel": r"$\Delta$ Si-Clu",
                    "suptitle": r"$\Delta$ Si-Clu:",
                    "title": rf"single $\mu$, {registry.get(det_mod).get_name(args.detname)}"
                    f" {' (no outliers)' if args.rm_outliers else ''}",
                }
                plot_track_param_hist_var_groups(
                    data,
                    threshold_outlier_detection,
                    args,
                    det_mod,
                    nominal_var_groups_labels,
                    group,
                )

    if args.detmods:
        for track_type in track_types:
            for var in ["D0", "Omega"]:
                nominal_diff_det_mods_labels = {
                    "suptitle": f"Diff DetMods: {var}",
                    "title": rf"single $\mu$, {track_type}"
                    f" {' (no outliers)' if args.rm_outliers else ''}",
                }
                plot_track_param_hist_diff_detmods(
                    data,
                    threshold_outlier_detection,
                    args,
                    f"{track_type}{var}",
                    nominal_diff_det_mods_labels,
                )


#############################################
# commands to access cov Matrix
#############################################
# import ROOT
# ROOT.gInterpreter.LoadFile("edm4hep/utils/cov_matrix_utils.h")
# # ...
# edm4hep.utils.get_cov_value(cov_m, edm4hep.TrackParams.d0, edm4hep.TrackParams.d0)
