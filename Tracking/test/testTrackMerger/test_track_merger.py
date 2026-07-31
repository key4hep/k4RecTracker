"""
Minimalistic test for the TrackMerger Gaudi processor.

What it does
------------
1. (run step) Creates a small EDM4hep input file containing:
   - InnerTracks : 2 tracks
   - OuterTracks : 2 tracks

   Track pair 0 is set up to MATCH  (|d0_diff| <= 0.5, |z0_diff| <= 2.5).
   Track pair 1 is set up to NOT match (large d0/z0 offsets).

   Then runs TrackMerger via test_track_merger_steer.py (k4run), which sets
   all 5 per-parameter tolerances (D0Tolerance, Z0Tolerance, PhiTolerance,
   OmegaTolerance, TanLambdaTolerance) to exercise their configurability.

2. (check step) Reads the output file and asserts that exactly 1 merged track
   was produced, that it carries hits from both parent tracks, and that it
   links back to both parent tracks.

The run and check steps are split into two subcommands so that CMake can run
them as separate tests with a fixture dependency between them, while sharing
all collection names/tolerances defined once in this file.

Usage
-----
    python3 test_track_merger.py run   --input input.root --output output.root
    python3 test_track_merger.py check --output output.root

Requirements
------------
  pip install podio edm4hep   (or load the key4hep stack)
  k4run must be on PATH
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


# ---------------------------------------------------------------------------
# Helper: build one TrackState with the minimum required fields
# ---------------------------------------------------------------------------
def make_track_state(location: int, d0: float, z0: float) -> edm4hep.TrackState:
    ts = edm4hep.TrackState()
    ts.location = location
    ts.D0 = d0
    ts.Z0 = z0
    ts.phi = 0.0
    ts.omega = 1e-4  # small curvature, arbitrary
    ts.tanLambda = 0.5
    return ts


# ---------------------------------------------------------------------------
# Helper: add a track with one TrackState and one dummy TrackerHit
# ---------------------------------------------------------------------------
def add_track(collection, hit_collection, location: int, d0: float, z0: float):
    hit = hit_collection.create()
    hit.setPosition(edm4hep.Vector3d(0.0, 0.0, 0.0))

    trk = collection.create()
    trk.addToTrackerHits(hit)
    trk.addToTrackStates(make_track_state(location, d0, z0))
    return trk


# ---------------------------------------------------------------------------
# Step 1 - write the input file
# ---------------------------------------------------------------------------
def write_input_file(path: str) -> None:
    writer = podio.root_io.Writer(path)

    frame = podio.Frame()

    inner_tracks = edm4hep.TrackCollection()
    outer_tracks = edm4hep.TrackCollection()

    inner_hits = edm4hep.TrackerHit3DCollection()
    outer_hits = edm4hep.TrackerHit3DCollection()

    # --- Pair 0: should MATCH ---
    # Inner track needs AtLastHit; Outer track needs AtFirstHit
    # Using identical d0/z0 -> diff = 0, well within tolerance
    d0_match, z0_match = 1.0, 5.0
    add_track(inner_tracks, inner_hits, edm4hep.TrackState.AtLastHit, d0_match, z0_match)
    add_track(outer_tracks, outer_hits, edm4hep.TrackState.AtFirstHit, d0_match, z0_match)

    # --- Pair 1: should NOT match ---
    # Offsets deliberately exceed both thresholds
    add_track(inner_tracks, inner_hits, edm4hep.TrackState.AtLastHit, 10.0, 100.0)
    add_track(outer_tracks, outer_hits, edm4hep.TrackState.AtFirstHit, 50.0, 500.0)

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
def check_output(output_file: str) -> None:
    reader = podio.root_io.Reader(output_file)

    frames = list(reader.get("events"))
    assert len(frames) == 1, f"Expected 1 event frame, got {len(frames)}"

    frame = frames[0]

    assert OUT_COLL in frame.getAvailableCollections(), (
        f"Output collection '{OUT_COLL}' not found in output file"
    )

    merged = frame.get(OUT_COLL)

    # --- core assertion: exactly one pair should have merged ---
    assert len(merged) == 1, (
        f"Expected 1 merged track (pair 0 matches, pair 1 does not), got {len(merged)}"
    )

    merged_track = merged[0]

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

    print("[check] All assertions passed.")
    print(f"        merged tracks : {len(merged)}")
    print(f"        parent links  : {len(parent_tracks)}")
    print(f"        merged hits   : {merged_hits}")


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

    args = parser.parse_args()

    if args.step == "setup":
        write_input_file(args.input)
    elif args.step == "check":
        check_output(args.output)


if __name__ == "__main__":
    main()
