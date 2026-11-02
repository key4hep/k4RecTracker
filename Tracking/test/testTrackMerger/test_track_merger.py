"""
Minimalistic test for the TrackMerger Gaudi processor.

What it does
------------
The test runs in one of two modes, selected with --mode:

absolute (default)
   - Pair 0 is set up to MATCH  (|d0_diff| <= 0.5, |z0_diff| <= 2.5).
   - Pair 1 is set up to NOT match (large d0/z0 offsets).

significance
   All track states carry covariances, and TrackMerger gets the magnetic
   field from field_only.xml (constant 2 T), so the pT significance is used.
   - Pair A: identical parameters                          -> MATCH
   - Pair B: d0 differs by 0.1 mm with sigma(d0) = 0.01 mm  -> NO match
     (would pass the absolute 0.5 mm cut, so this proves significances are used)
   - Pair C: same |omega| but opposite sign (charge)       -> NO match
   - Pair D: phi = +3.13 vs -3.13 (across +-pi)            -> MATCH
   The pairs are well separated in d0/z0 so that they cannot cross-match.

1. (setup step) Creates a small EDM4hep input file with InnerTracks and
   OuterTracks for the selected mode.

   TrackMerger itself is run separately by CTest via
   `k4run test_track_merger_steer.py --mode <mode>`, which sets all 5
   per-parameter cuts of that mode to exercise their configurability.

2. (check step) Reads the output file and asserts that exactly the expected
   pairs were merged, that each merged track carries hits from both parent
   tracks, and that it links back to both parent tracks.

The setup and check steps are split into two subcommands so that CMake can
run them (and the k4run step in between) as separate tests chained via
fixtures, while sharing all collection names defined once in this file.

Usage
-----
    python3 test_track_merger.py setup --mode <mode> --input input.root
    k4run test_track_merger_steer.py   --mode <mode> --input input.root --output output.root
    python3 test_track_merger.py check --mode <mode> --output output.root

Requirements
------------
  pip install podio edm4hep   (or load the key4hep stack)
"""

import argparse

import edm4hep
import podio

# ---------------------------------------------------------------------------
# Constants that mirror the TrackMerger defaults.
# Must be kept in sync with the constants in test_track_merger_steer.py.
# ---------------------------------------------------------------------------
INNER_COLL = "InnerTracks"
OUTER_COLL = "OuterTracks"
OUT_COLL = "MyCandidateMergedTracks"


# Default variances used in significance mode (sigma = 0.1 mm, 0.01 rad, ...).
DEFAULT_VARIANCES = {
    edm4hep.TrackParams.d0: 1e-2,
    edm4hep.TrackParams.z0: 1e-2,
    edm4hep.TrackParams.phi: 1e-4,
    edm4hep.TrackParams.omega: 1e-10,
    edm4hep.TrackParams.tanLambda: 1e-4,
}


# ---------------------------------------------------------------------------
# Helper: build one TrackState with the minimum required fields
# ---------------------------------------------------------------------------
def make_track_state(
    location: int,
    d0: float,
    z0: float,
    phi: float = 0.0,
    omega: float = 1e-4,
    variances: dict | None = None,
) -> edm4hep.TrackState:
    ts = edm4hep.TrackState()
    ts.location = location
    ts.D0 = d0
    ts.Z0 = z0
    ts.phi = phi
    ts.omega = omega
    ts.tanLambda = 0.5
    for param, var in (variances or {}).items():
        ts.covMatrix.setValue(var, param, param)
    return ts


# ---------------------------------------------------------------------------
# Helper: add a track with one TrackState and one dummy TrackerHit
# ---------------------------------------------------------------------------
def add_track(collection, hit_collection, location: int, d0: float, z0: float, **state_kwargs):
    hit = hit_collection.create()
    hit.setPosition(edm4hep.Vector3d(0.0, 0.0, 0.0))

    trk = collection.create()
    trk.addToTrackerHits(hit)
    trk.addToTrackStates(make_track_state(location, d0, z0, **state_kwargs))
    return trk


def add_pair(inner_tracks, inner_hits, outer_tracks, outer_hits, inner_kwargs, outer_kwargs):
    # Inner track needs AtLastHit; Outer track needs AtFirstHit
    add_track(inner_tracks, inner_hits, edm4hep.TrackState.AtLastHit, **inner_kwargs)
    add_track(outer_tracks, outer_hits, edm4hep.TrackState.AtFirstHit, **outer_kwargs)


def fill_absolute(inner_tracks, inner_hits, outer_tracks, outer_hits):
    tracks = (inner_tracks, inner_hits, outer_tracks, outer_hits)
    # --- Pair 0: should MATCH ---
    # Using identical d0/z0 -> diff = 0, well within tolerance
    add_pair(*tracks, dict(d0=1.0, z0=5.0), dict(d0=1.0, z0=5.0))
    # --- Pair 1: should NOT match ---
    # Offsets deliberately exceed both thresholds
    add_pair(*tracks, dict(d0=10.0, z0=100.0), dict(d0=50.0, z0=500.0))


def fill_significance(inner_tracks, inner_hits, outer_tracks, outer_hits):
    tracks = (inner_tracks, inner_hits, outer_tracks, outer_hits)
    common = dict(omega=1e-3, variances=DEFAULT_VARIANCES)
    # --- Pair A: identical -> MATCH ---
    add_pair(*tracks, dict(d0=1.0, z0=0.0, **common), dict(d0=1.0, z0=0.0, **common))
    # --- Pair B: 0.1 mm in d0 with sigma(d0) = 0.01 mm -> ~7 sigma -> NO match ---
    precise_d0 = dict(DEFAULT_VARIANCES)
    precise_d0[edm4hep.TrackParams.d0] = 1e-4
    add_pair(
        *tracks,
        dict(d0=20.0, z0=100.0, omega=1e-3, variances=precise_d0),
        dict(d0=20.1, z0=100.0, omega=1e-3, variances=precise_d0),
    )
    # --- Pair C: opposite charge, same curvature -> NO match ---
    add_pair(
        *tracks,
        dict(d0=40.0, z0=200.0, omega=1e-3, variances=DEFAULT_VARIANCES),
        dict(d0=40.0, z0=200.0, omega=-1e-3, variances=DEFAULT_VARIANCES),
    )
    # --- Pair D: across phi = +-pi, wrapped diff ~0.023 rad -> ~1.6 sigma -> MATCH ---
    add_pair(
        *tracks,
        dict(d0=60.0, z0=300.0, phi=3.13, **common),
        dict(d0=60.0, z0=300.0, phi=-3.13, **common),
    )


# d0 of the inner tracks of the pairs that are expected to merge, per mode
EXPECTED_MERGED_D0 = {
    "absolute": [1.0],
    "significance": [1.0, 60.0],
}


# ---------------------------------------------------------------------------
# Step 1 - write the input file
# ---------------------------------------------------------------------------
def write_input_file(path: str, mode: str) -> None:
    writer = podio.root_io.Writer(path)

    frame = podio.Frame()

    inner_tracks = edm4hep.TrackCollection()
    outer_tracks = edm4hep.TrackCollection()

    inner_hits = edm4hep.TrackerHit3DCollection()
    outer_hits = edm4hep.TrackerHit3DCollection()

    fill = fill_absolute if mode == "absolute" else fill_significance
    fill(inner_tracks, inner_hits, outer_tracks, outer_hits)

    frame.put(inner_tracks, INNER_COLL)
    frame.put(outer_tracks, OUTER_COLL)
    frame.put(inner_hits, "InnerHits")
    frame.put(outer_hits, "OuterHits")

    writer.write_frame(frame, "events")
    writer.finish()
    print(f"[setup] Wrote input file: {path}")


# ---------------------------------------------------------------------------
# Step 2 - read back & assert
# ---------------------------------------------------------------------------
def check_output(output_file: str, mode: str) -> None:
    reader = podio.root_io.Reader(output_file)

    frames = list(reader.get("events"))
    assert len(frames) == 1, f"Expected 1 event frame, got {len(frames)}"

    frame = frames[0]

    assert OUT_COLL in frame.getAvailableCollections(), (
        f"Output collection '{OUT_COLL}' not found in output file"
    )

    merged = frame.get(OUT_COLL)

    # --- core assertion: exactly the expected pairs should have merged ---
    expected_d0 = EXPECTED_MERGED_D0[mode]
    assert len(merged) == len(expected_d0), (
        f"[{mode}] Expected {len(expected_d0)} merged track(s), got {len(merged)}"
    )

    merged_d0 = []
    for merged_track in merged:
        # --- the merged track must link back to both parent tracks ---
        parent_tracks = list(merged_track.getTracks())
        assert len(parent_tracks) == 2, (
            f"Merged track should carry 2 sub-tracks (inner + outer), got {len(parent_tracks)}"
        )

        # --- it must aggregate hits from both parents ---
        total_parent_hits = sum(len(list(p.getTrackerHits())) for p in parent_tracks)
        merged_hits = len(list(merged_track.getTrackerHits()))
        assert merged_hits == total_parent_hits, (
            f"Merged track should have {total_parent_hits} hits (1 from inner + 1 from outer), "
            f"got {merged_hits}"
        )

        merged_d0.append(round(parent_tracks[0].getTrackStates(0).D0, 3))

    # --- the right pairs must have merged ---
    assert sorted(merged_d0) == sorted(expected_d0), (
        f"[{mode}] Expected merged pairs with d0 {sorted(expected_d0)}, got {sorted(merged_d0)}"
    )

    print(f"[check] All assertions passed ({mode} mode).")
    print(f"        merged tracks : {len(merged)}")
    print(f"        merged pairs  : d0 = {sorted(merged_d0)}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="step", required=True)

    setup_parser = subparsers.add_parser("setup", help="Write the input file")
    setup_parser.add_argument("--input", default="input.root")

    check_parser = subparsers.add_parser("check", help="Check the TrackMerger output file")
    check_parser.add_argument("--output", default="output.root")

    for sub in (setup_parser, check_parser):
        sub.add_argument("--mode", choices=list(EXPECTED_MERGED_D0), default="absolute")

    args = parser.parse_args()

    if args.step == "setup":
        write_input_file(args.input, args.mode)
    elif args.step == "check":
        check_output(args.output, args.mode)


if __name__ == "__main__":
    main()
