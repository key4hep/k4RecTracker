"""
Minimalistic test for the SimTrackerHitCellMerger Gaudi transformer.

What it does
------------
1. (setup step) Creates a small EDM4hep input file containing three MCParticles and one
   SimTrackerHitCollection with hits in two cells:

     cell 1 : particle 0 (the "primary", lowest Geant4 track number) with 2 steps of 1.0 and 2.0 mm
              particle 1                                             with 1 step  of 0.5 mm
              particle 2                                             with 1 step  of 0.25 mm
     cell 2 : particle 1                                             with 1 step  of 4.0 mm

   The hits are deliberately written in a scrambled order so that the per-cell and per-track
   grouping, the cellID ordering of the output and the "earliest hit is the representative one"
   rule are all exercised.

   SimTrackerHitCellMerger itself is run separately by CTest via
   `k4run test_simtrackerhit_cell_merger_steer.py`, which runs one instance per value of the
   MultipleTrackHandling property so that all three are covered in a single job.

2. (check step) Reads the output file and asserts the accumulated path lengths, energy deposits,
   MCParticle relations and representative-hit times of each of the three output collections.

Usage
-----
    python3 test_simtrackerhit_cell_merger.py setup   --input input.root
    k4run test_simtrackerhit_cell_merger_steer.py     --input input.root --output output.root
    python3 test_simtrackerhit_cell_merger.py check   --output output.root

Requirements
------------
  pip install podio edm4hep   (or load the key4hep stack)
"""

import argparse

import edm4hep
import podio

# ---------------------------------------------------------------------------
# Constants that mirror the steering file.
# Must be kept in sync with the constants in test_simtrackerhit_cell_merger_steer.py.
# ---------------------------------------------------------------------------
PARTICLE_COLL = "MCParticles"
INPUT_COLL = "SimTrackerHits"
OUT_COLL_ALL = "MergedHitsAll"
OUT_COLL_PRIMARY = "MergedHitsPrimaryOnly"
OUT_COLL_PER_TRACK = "MergedHitsPerTrack"

CELL_A = 1
CELL_B = 2

# (cellID, particle index, pathLength [mm], eDep [GeV], time [ns]), in scrambled order
HITS = [
    (CELL_A, 0, 2.0, 0.002, 2.0),
    (CELL_B, 1, 4.0, 0.004, 3.0),
    (CELL_A, 2, 0.25, 0.00025, 0.5),
    (CELL_A, 0, 1.0, 0.001, 1.0),
    (CELL_A, 1, 0.5, 0.0005, 1.5),
]


# ---------------------------------------------------------------------------
# Step 1 - write the input file
# ---------------------------------------------------------------------------
def write_input_file(path: str) -> None:
    writer = podio.root_io.Writer(path)

    frame = podio.Frame()

    particles = edm4hep.MCParticleCollection()
    # DD4hep fills the MCParticle collection ordered by increasing Geant4 trackID, so the index
    # within the collection is what SimTrackerHitCellMerger uses to rank the contributors.
    # Particle 0 plays the role of the primary muon, 1 and 2 those of two later secondaries.
    for pdg in (13, 11, 11):
        particle = particles.create()
        particle.setPDG(pdg)
        particle.setMomentum(edm4hep.Vector3d(1.0, 0.0, 0.0))

    hits = edm4hep.SimTrackerHitCollection()
    for cell_id, particle_index, path_length, edep, time in HITS:
        hit = hits.create()
        hit.setCellID(cell_id)
        hit.setPathLength(path_length)
        hit.setEDep(edep)
        hit.setTime(time)
        hit.setPosition(edm4hep.Vector3d(time, 0.0, 0.0))
        hit.setMomentum(edm4hep.Vector3f(1.0, 0.0, 0.0))
        hit.setParticle(particles[particle_index])

    frame.put(particles, PARTICLE_COLL)
    frame.put(hits, INPUT_COLL)

    writer.write_frame(frame, "events")
    writer.finish()
    print(f"[setup] Wrote input file: {path}")


# ---------------------------------------------------------------------------
# Helper: compare one merged hit against its expectation
# ---------------------------------------------------------------------------
def check_hit(
    coll_name: str,
    index: int,
    hit,
    cell_id: int,
    path_length: float,
    edep: float,
    particle_index: int,
    time: float,
) -> None:
    prefix = f"{coll_name}[{index}]"
    assert hit.getCellID() == cell_id, (
        f"{prefix}: expected cellID {cell_id}, got {hit.getCellID()}"
    )
    assert abs(hit.getPathLength() - path_length) < 1e-5, (
        f"{prefix}: expected path length {path_length} mm, got {hit.getPathLength()}"
    )
    assert abs(hit.getEDep() - edep) < 1e-8, (
        f"{prefix}: expected eDep {edep} GeV, got {hit.getEDep()}"
    )
    assert hit.getParticle().getObjectID().index == particle_index, (
        f"{prefix}: expected MCParticle {particle_index}, got {hit.getParticle().getObjectID().index}"
    )
    # The representative hit is the earliest one that was summed, and it provides the time
    assert abs(hit.getTime() - time) < 1e-5, (
        f"{prefix}: expected time {time} ns, got {hit.getTime()}"
    )


# ---------------------------------------------------------------------------
# Step 2 - read back & assert
# ---------------------------------------------------------------------------
def check_output(output_file: str) -> None:
    reader = podio.root_io.Reader(output_file)

    frames = list(reader.get("events"))
    assert len(frames) == 1, f"Expected 1 event frame, got {len(frames)}"

    frame = frames[0]
    available = frame.getAvailableCollections()
    for coll_name in (OUT_COLL_ALL, OUT_COLL_PRIMARY, OUT_COLL_PER_TRACK):
        assert coll_name in available, f"Output collection '{coll_name}' not found in output file"

    # --- MultipleTrackHandling = "All" -------------------------------------
    # One hit per cell, summing every track. Cell 1 gets 1.0 + 2.0 + 0.5 + 0.25 mm, is attributed to
    # the most primary contributor (particle 0) and takes its time from the globally earliest hit of
    # the cell, which belongs to particle 2.
    merged_all = frame.get(OUT_COLL_ALL)
    assert len(merged_all) == 2, (
        f"'All' should give one hit per cell, i.e. 2, got {len(merged_all)}"
    )
    check_hit(OUT_COLL_ALL, 0, merged_all[0], CELL_A, 3.75, 0.00375, 0, 0.5)
    check_hit(OUT_COLL_ALL, 1, merged_all[1], CELL_B, 4.0, 0.004, 1, 3.0)

    # --- MultipleTrackHandling = "PrimaryOnly" -----------------------------
    # One hit per cell, but only the steps of the lowest-track-number contributor are summed, so
    # cell 1 keeps only particle 0's 1.0 + 2.0 mm and the 0.5 + 0.25 mm of the secondaries are dropped.
    merged_primary = frame.get(OUT_COLL_PRIMARY)
    assert len(merged_primary) == 2, (
        f"'PrimaryOnly' should give one hit per cell, i.e. 2, got {len(merged_primary)}"
    )
    check_hit(OUT_COLL_PRIMARY, 0, merged_primary[0], CELL_A, 3.0, 0.003, 0, 1.0)
    check_hit(OUT_COLL_PRIMARY, 1, merged_primary[1], CELL_B, 4.0, 0.004, 1, 3.0)

    # --- MultipleTrackHandling = "PerTrack" --------------------------------
    # Same data type, but cell 1 now appears three times, once per contributing track and ordered by
    # increasing Geant4 track number.
    merged_per_track = frame.get(OUT_COLL_PER_TRACK)
    assert len(merged_per_track) == 4, (
        f"'PerTrack' should give 3 hits for cell {CELL_A} and 1 for cell {CELL_B}, i.e. 4, "
        f"got {len(merged_per_track)}"
    )
    check_hit(OUT_COLL_PER_TRACK, 0, merged_per_track[0], CELL_A, 3.0, 0.003, 0, 1.0)
    check_hit(OUT_COLL_PER_TRACK, 1, merged_per_track[1], CELL_A, 0.5, 0.0005, 1, 1.5)
    check_hit(OUT_COLL_PER_TRACK, 2, merged_per_track[2], CELL_A, 0.25, 0.00025, 2, 0.5)
    check_hit(OUT_COLL_PER_TRACK, 3, merged_per_track[3], CELL_B, 4.0, 0.004, 1, 3.0)

    print("[check] All assertions passed.")
    print(f"        {OUT_COLL_ALL}       : {len(merged_all)} hits")
    print(f"        {OUT_COLL_PRIMARY}   : {len(merged_primary)} hits")
    print(f"        {OUT_COLL_PER_TRACK} : {len(merged_per_track)} hits")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="step", required=True)

    setup_parser = subparsers.add_parser("setup", help="Write the input file")
    setup_parser.add_argument("--input", default="input.root")

    check_parser = subparsers.add_parser(
        "check", help="Check the SimTrackerHitCellMerger output file"
    )
    check_parser.add_argument("--output", default="output.root")

    args = parser.parse_args()

    if args.step == "setup":
        write_input_file(args.input)
    elif args.step == "check":
        check_output(args.output)


if __name__ == "__main__":
    main()
