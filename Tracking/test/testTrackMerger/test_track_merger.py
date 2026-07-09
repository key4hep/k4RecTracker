"""
Minimalistic test for the TrackMerger Gaudi processor.

What it does
------------
1. Creates a small EDM4hep input file containing:
   - SiTracksCT     : 2 tracks
   - ClupatraTracks : 2 tracks

   Track pair 0 is set up to MATCH  (|d0_diff| <= 0.5, |z0_diff| <= 2.5).
   Track pair 1 is set up to NOT match (large d0/z0 offsets).

2. Runs TrackMerger via a minimal Gaudi job (k4run).

3. Reads the output file and asserts that exactly 1 merged track was produced,
   that it carries hits from both parent tracks, and that it links back to both
   parent tracks.

Requirements
------------
  pip install podio edm4hep   (or load the key4hep stack)
  k4run must be on PATH
"""

import subprocess
import tempfile
from pathlib import Path

import edm4hep
import podio
from edm4hep import TrackState as TS

# ---------------------------------------------------------------------------
# Constants that mirror the TrackMerger defaults
# ---------------------------------------------------------------------------
SI_COLL = "SiTracksCT"
CLU_COLL = "ClupatraTracks"
OUT_COLL = "MyCandidateMergedTracks"

# Match thresholds (from TrackMerger.cpp)
D0_TOL = 0.5
Z0_TOL = 2.5


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

    si_tracks = edm4hep.TrackCollection()
    clu_tracks = edm4hep.TrackCollection()

    si_hits = edm4hep.TrackerHit3DCollection()
    clu_hits = edm4hep.TrackerHit3DCollection()

    # --- Pair 0: should MATCH ---
    # SiTrack needs AtLastHit; CluTrack needs AtFirstHit
    # Using identical d0/z0 -> diff = 0, well within tolerance
    d0_match, z0_match = 1.0, 5.0
    add_track(si_tracks, si_hits, TS.AtLastHit, d0_match, z0_match)
    add_track(clu_tracks, clu_hits, TS.AtFirstHit, d0_match, z0_match)

    # --- Pair 1: should NOT match ---
    # Offsets deliberately exceed both thresholds
    add_track(si_tracks, si_hits, TS.AtLastHit, 10.0, 100.0)
    add_track(clu_tracks, clu_hits, TS.AtFirstHit, 50.0, 500.0)

    frame.put(si_tracks, SI_COLL)
    frame.put(clu_tracks, CLU_COLL)
    frame.put(si_hits, "SiHits")
    frame.put(clu_hits, "CluHits")

    writer.write_frame(frame, "events")
    writer.finish()
    print(f"[setup] Wrote input file: {path}")


# ---------------------------------------------------------------------------
# Step 2 - write a minimal Gaudi steering file on-the-fly
# ---------------------------------------------------------------------------
STEERING_TEMPLATE = """\
from Configurables import TrackMerger
from k4FWCore import ApplicationMgr, IOSvc

iosvc = IOSvc()
iosvc.Input  = "{input}"
iosvc.Output = "{output}"

merger = TrackMerger("TrackMerger",
    InputInnerTracks  = "{si_coll}",
    InputOuterTracks = "{clu_coll}",
    OutTracks      = "{out_coll}",
    Greedy         = True,
)

ApplicationMgr(
    TopAlg  = [merger],
    EvtSel  = "NONE",
    EvtMax  = 1,
    ExtSvc  = [iosvc],
)
"""


def write_steering_file(path: str, input_file: str, output_file: str) -> None:
    content = STEERING_TEMPLATE.format(
        input=input_file,
        output=output_file,
        si_coll=SI_COLL,
        clu_coll=CLU_COLL,
        out_coll=OUT_COLL,
    )
    Path(path).write_text(content)
    print(f"[setup] Wrote steering file: {path}")


# ---------------------------------------------------------------------------
# Step 3 - run Gaudi
# ---------------------------------------------------------------------------
def run_gaudi(steering_file: str) -> None:
    result = subprocess.run(
        ["k4run", steering_file],
        capture_output=True,
        text=True,
    )
    print(result.stdout[-3000:])  # tail to keep output readable
    if result.returncode != 0:
        print(result.stderr[-2000:])
        raise RuntimeError(f"k4run failed (exit {result.returncode})")
    print("[run] Gaudi job finished successfully.")


# ---------------------------------------------------------------------------
# Step 4 - read back & assert
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
        f"Merged track should carry 2 sub-tracks (SiTrack + CluTrack), got {len(parent_tracks)}"
    )

    # --- it must aggregate hits from both parents ---
    total_parent_hits = sum(len(list(p.getTrackerHits())) for p in parent_tracks)
    merged_hits = len(list(merged_track.getTrackerHits()))
    assert merged_hits == total_parent_hits, (
        f"Merged track should have {total_parent_hits} hits (1 from Si + 1 from Clu), "
        f"got {merged_hits}"
    )

    print("[check] All assertions passed.")
    print(f"        merged tracks : {len(merged)}")
    print(f"        parent links  : {len(parent_tracks)}")
    print(f"        merged hits   : {merged_hits}")


# ---------------------------------------------------------------------------
# Entry point (also works as a standalone script)
# ---------------------------------------------------------------------------
def test_track_merger():
    with tempfile.TemporaryDirectory() as tmpdir:
        input_file = str(Path(tmpdir) / "input.root")
        output_file = str(Path(tmpdir) / "output.root")
        steering_file = str(Path(tmpdir) / "steer.py")

        write_input_file(input_file)
        write_steering_file(steering_file, input_file, output_file)
        run_gaudi(steering_file)
        check_output(output_file)


if __name__ == "__main__":
    test_track_merger()
