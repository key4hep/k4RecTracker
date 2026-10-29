#!/usr/bin/env python3
# Produces the comparison plots for the two tracking pipelines (GGTF+Genfit vs
# ConformalTracking+RefitFinal) on the IDEA_o1_v04 samples.
# Reuses the analysis/matching from compare_tracking.py; files must already
# have been produced by runIDEA_GGTF_Genfit.py and runIDEA_ConformalTracking_Refit.py.
#
# Usage: python3 plot_tracking.py [fileA.root] [fileB.root] [--truth gun|stable]
# Saves: tracking_residuals.png, tracking_hits_chi2.png, tracking_efficiency.png

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import compare_tracking
from compare_tracking import analyze

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("A", nargs="?", default="tracks_ggtf_genfit.edm4hep.root")
ap.add_argument("B", nargs="?", default="tracks_ct_refit.edm4hep.root")
ap.add_argument("--truth", choices=["gun", "stable"], default="gun")
ap.add_argument("--prefix", default="tracking", help="output file prefix")
args = ap.parse_args()

compare_tracking.TRUTH_MODE = args.truth

outA = analyze(args.A, "FittedTracks")
outB = analyze(args.B, "SiTracksCT_Refitted")

LABEL_A = "GGTF + Genfit"
LABEL_B = "ConformalTracking + Refit"
COL_A, COL_B = "tab:orange", "tab:blue"


def overlay(ax, va, vb, bins, xlabel, tag):
    ax.hist(
        va,
        bins=bins,
        histtype="stepfilled",
        alpha=0.45,
        color=COL_A,
        label=f"{LABEL_A} (n={len(va)})",
    )
    ax.hist(
        vb,
        bins=bins,
        histtype="stepfilled",
        alpha=0.45,
        color=COL_B,
        label=f"{LABEL_B} (n={len(vb)})",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("matched tracks")
    ax.legend(fontsize=8)
    ax.set_title(tag, fontsize=10)


def rms(v):
    return float(np.std(v, ddof=1)) if len(v) > 1 else 0.0


# ---------------- residuals ----------------
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
raw_a, raw_b = outA["raw"], outB["raw"]
overlay(
    axes[0, 0],
    raw_a["rel_pt"],
    raw_b["rel_pt"],
    np.linspace(-0.08, 0.08, 41),
    r"$(p_T^{\rm reco}-p_T^{\rm MC})/p_T^{\rm MC}$",
    "transverse momentum",
)
overlay(
    axes[0, 1],
    [x * 1e3 for x in raw_a["phi"]],
    [x * 1e3 for x in raw_b["phi"]],
    np.linspace(-10, 10, 41),
    r"$\Delta\phi$ [mrad]",
    "azimuthal angle",
)
overlay(
    axes[0, 2],
    [x * 1e3 for x in raw_a["tanLambda"]],
    [x * 1e3 for x in raw_b["tanLambda"]],
    np.linspace(-10, 10, 41),
    r"$\Delta\tan\lambda$ [$10^{-3}$]",
    "track slope",
)
overlay(
    axes[1, 0],
    raw_a["d0"],
    raw_b["d0"],
    np.linspace(-0.1, 0.1, 41),
    r"$d_0^{\rm reco}-d_0^{\rm MC}$ [mm]",
    "transverse impact parameter",
)
overlay(
    axes[1, 1],
    raw_a["z0"],
    raw_b["z0"],
    np.linspace(-0.2, 0.2, 41),
    r"$z_0^{\rm reco}-z_0^{\rm MC}$ [mm]",
    "longitudinal impact parameter",
)
ax = axes[1, 2]
names = ["rel pT [%]", r"$\phi$ [mrad]", r"$\tan\lambda$ [$10^{-3}$]", "d0 [mm]", "z0 [mm]"]
rms_a = [
    rms(raw_a["rel_pt"]) * 100,
    rms(raw_a["phi"]) * 1e3,
    rms(raw_a["tanLambda"]) * 1e3,
    rms(raw_a["d0"]),
    rms(raw_a["z0"]),
]
rms_b = [
    rms(raw_b["rel_pt"]) * 100,
    rms(raw_b["phi"]) * 1e3,
    rms(raw_b["tanLambda"]) * 1e3,
    rms(raw_b["d0"]),
    rms(raw_b["z0"]),
]
x = np.arange(len(names))
ax.bar(x - 0.2, rms_a, width=0.4, color=COL_A, label=LABEL_A)
ax.bar(x + 0.2, rms_b, width=0.4, color=COL_B, label=LABEL_B)
ax.set_yscale("log")
ax.set_xticks(x, names, fontsize=8)
ax.set_ylabel("RMS")
ax.set_title("resolution summary", fontsize=10)
ax.legend(fontsize=8)
fig.suptitle(
    f"IDEA_o1_v04, truth mode {args.truth} — residual to generator truth (matched tracks)",
    fontsize=12,
)
fig.tight_layout()
fig.savefig(f"{args.prefix}_residuals.png", dpi=150)

# ---------------- hits and chi2 ----------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
max_hits = max(raw_a["hits"] + raw_b["hits"])
overlay(
    axes[0],
    raw_a["hits"],
    raw_b["hits"],
    np.arange(0, max_hits + 2) - 0.5,
    "hits per track",
    "track hit content",
)
overlay(
    axes[1],
    raw_a["chi2ndf"],
    raw_b["chi2ndf"],
    np.linspace(0, 8, 41),
    r"$\chi^2 / ndf$",
    "fit quality",
)
ax = axes[2]
colls = list(outA["hit_fraction"].keys())
fa = [outA["hit_fraction"][c] for c in colls]
fb = [outB["hit_fraction"][c] for c in colls]
x = np.arange(len(colls))
ax.bar(x - 0.2, fa, width=0.4, color=COL_A, label=LABEL_A)
ax.bar(x + 0.2, fb, width=0.4, color=COL_B, label=LABEL_B)
ax.set_xticks(x, [c.replace("Digis", "") for c in colls], fontsize=8)
ax.set_ylabel("fraction of attached planar hits")
ax.set_title("hit provenance on tracks", fontsize=10)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{args.prefix}_hits_chi2.png", dpi=150)

# ---------------- per-event efficiency ----------------
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
for ax, out, label, col in ((axes[0], outA, LABEL_A, COL_A), (axes[1], outB, LABEL_B, COL_B)):
    evs = out["events"]
    idx = [e["event"] for e in evs]
    eff = [e["matched"] / e["muons"] if e["muons"] else 0 for e in evs]
    mu = [e["muons"] for e in evs]
    trk = [e["tracks"] for e in evs]
    ax.bar(idx, eff, color=col, alpha=0.7)
    ax.set_xticks(idx)
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("event")
    ax.set_ylabel("matched / truth particles")
    tot = out["counters"]
    ax.set_title(
        f"{label}\noverall: {tot['matched']}/{tot['muons']} = "
        f"{tot['matched'] / tot['muons']:.2f}, tracks/ev={tot['tracks'] / tot['events']:.1f} "
        f"(+{tot['failed']} failed)",
        fontsize=9,
    )
    ax2 = ax.twinx()
    ax2.plot(idx, trk, "k.--", lw=1, label="tracks (incl. failed)")
    ax2.plot(idx, mu, "rs:", lw=1, label="truth particles")
    ax2.set_ylabel("count")
    ax2.set_ylim(0, max(max(trk), max(mu)) + 2)
    ax2.legend(fontsize=8, loc="lower right")
fig.suptitle("truth matching efficiency per event", fontsize=12)
fig.tight_layout()
fig.savefig(f"{args.prefix}_efficiency.png", dpi=150)

print(
    f"saved: {args.prefix}_residuals.png {args.prefix}_hits_chi2.png {args.prefix}_efficiency.png"
)
