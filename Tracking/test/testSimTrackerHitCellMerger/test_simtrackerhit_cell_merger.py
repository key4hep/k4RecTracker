"""
Minimalistic test for the SimTrackerHitCellMerger Gaudi multi transformer.

What it does
------------
1. (setup step) Creates a small EDM4hep input file containing three MCParticles and one
   SimTrackerHitCollection with hits in two cells:

     cell 1 : particle 0 (the most primary one) with 2 steps of 1.0 and 2.0 mm
              particle 1                        with 1 step  of 0.5 mm
              particle 2                        with 1 step  of 0.25 mm
     cell 2 : particle 1                        with 1 step  of 4.0 mm

   So cell 1 is crossed by three tracks and cell 2 by a single one. The hits are deliberately
   written in a scrambled order so that the per-cell and per-track grouping, the cellID ordering of
   the output and the "earliest hit is the representative one" rule are all exercised. Each hit
   carries a position and a momentum derived from its time, so that averaging them gives a
   predictable answer.

   SimTrackerHitCellMerger itself is run separately by CTest via
   `k4run test_simtrackerhit_cell_merger_steer.py`, which runs one instance per value of the
   MultipleTrackHandling property, plus one more for the "Average" RepresentativeKinematics, so
   that every choice of both properties is covered in a single job.

2. (check step) Reads the output file and asserts the accumulated path lengths, energy deposits,
   MCParticle relations and representative kinematics of each output collection.

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
OUT_COLL_SINGLE_TRACK = "MergedHitsSingleTrackCells"
OUT_COLL_AVERAGE = "MergedHitsAllEarliestHit"

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


def hit_position(time: float) -> tuple:
    """Position a hit at a point that is a unique function of its time."""
    return (time, 2.0 * time, 3.0 * time)


def hit_momentum(time: float) -> tuple:
    """Give a hit a momentum that is a unique function of its time."""
    return (time, 0.0, 0.0)


# ---------------------------------------------------------------------------
# Step 1 - write the input file
# ---------------------------------------------------------------------------
def write_input_file(path: str) -> None:
    writer = podio.root_io.Writer(path)

    frame = podio.Frame()

    particles = edm4hep.MCParticleCollection()
    # None of these is flagged as created in simulation, so SimTrackerHitCellMerger ranks them by
    # their index within the collection, which DD4hep writes in order of increasing Geant4 trackID.
    # Particle 0 therefore plays the role of the primary muon, 1 and 2 those of two later tracks.
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
        hit.setPosition(edm4hep.Vector3d(*hit_position(time)))
        hit.setMomentum(edm4hep.Vector3f(*hit_momentum(time)))
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
    """Assert the accumulated quantities, the MCParticle relation and the representative kinematics.

    `time` is the expected time of the merged hit; its position and momentum are expected to be the
    ones that go with that time, which holds both for a copy of a single input hit and for an
    average, because both are linear in the time.
    """
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
    assert abs(hit.getTime() - time) < 1e-5, (
        f"{prefix}: expected time {time} ns, got {hit.getTime()}"
    )

    expected_position = hit_position(time)
    position = (hit.getPosition().x, hit.getPosition().y, hit.getPosition().z)
    assert all(abs(a - b) < 1e-4 for a, b in zip(position, expected_position)), (
        f"{prefix}: expected position {expected_position} mm, got {position}"
    )

    expected_momentum = hit_momentum(time)
    momentum = (hit.getMomentum().x, hit.getMomentum().y, hit.getMomentum().z)
    assert all(abs(a - b) < 1e-4 for a, b in zip(momentum, expected_momentum)), (
        f"{prefix}: expected momentum {expected_momentum} GeV, got {momentum}"
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
    for coll_name in (
        OUT_COLL_ALL,
        OUT_COLL_PRIMARY,
        OUT_COLL_PER_TRACK,
        OUT_COLL_SINGLE_TRACK,
        OUT_COLL_AVERAGE,
    ):
        assert coll_name in available, f"Output collection '{coll_name}' not found in output file"

    # --- MultipleTrackHandling = "SumAll" ----------------------------------
    # One hit per cell, summing every track. Cell 1 gets 1.0 + 2.0 + 0.5 + 0.25 mm, is attributed to
    # the most primary contributor (particle 0) and takes its kinematics from the globally earliest
    # hit of the cell, which belongs to particle 2.
    merged_all = frame.get(OUT_COLL_ALL)
    assert len(merged_all) == 2, (
        f"'SumAll' should give one hit per cell, i.e. 2, got {len(merged_all)}"
    )
    check_hit(OUT_COLL_ALL, 0, merged_all[0], CELL_A, 3.75, 0.00375, 0, 0.5)
    check_hit(OUT_COLL_ALL, 1, merged_all[1], CELL_B, 4.0, 0.004, 1, 3.0)

    # --- MultipleTrackHandling = "PrimaryOnly" -----------------------------
    # One hit per cell, but only the steps of the most primary contributor are summed, so cell 1
    # keeps only particle 0's 1.0 + 2.0 mm and the 0.5 + 0.25 mm of the other two tracks are dropped.
    merged_primary = frame.get(OUT_COLL_PRIMARY)
    assert len(merged_primary) == 2, (
        f"'PrimaryOnly' should give one hit per cell, i.e. 2, got {len(merged_primary)}"
    )
    check_hit(OUT_COLL_PRIMARY, 0, merged_primary[0], CELL_A, 3.0, 0.003, 0, 1.0)
    check_hit(OUT_COLL_PRIMARY, 1, merged_primary[1], CELL_B, 4.0, 0.004, 1, 3.0)

    # --- MultipleTrackHandling = "PerTrack" --------------------------------
    # Same data type, but cell 1 now appears three times, once per contributing track and ordered
    # from the most to the least primary one.
    merged_per_track = frame.get(OUT_COLL_PER_TRACK)
    assert len(merged_per_track) == 4, (
        f"'PerTrack' should give 3 hits for cell {CELL_A} and 1 for cell {CELL_B}, i.e. 4, "
        f"got {len(merged_per_track)}"
    )
    check_hit(OUT_COLL_PER_TRACK, 0, merged_per_track[0], CELL_A, 3.0, 0.003, 0, 1.0)
    check_hit(OUT_COLL_PER_TRACK, 1, merged_per_track[1], CELL_A, 0.5, 0.0005, 1, 1.5)
    check_hit(OUT_COLL_PER_TRACK, 2, merged_per_track[2], CELL_A, 0.25, 0.00025, 2, 0.5)
    check_hit(OUT_COLL_PER_TRACK, 3, merged_per_track[3], CELL_B, 4.0, 0.004, 1, 3.0)

    # --- MultipleTrackHandling = "SkipMultiTrackCells" ---------------------
    # Cell 1 was crossed by three tracks and is dropped entirely; only the unambiguous cell 2 survives.
    merged_single_track = frame.get(OUT_COLL_SINGLE_TRACK)
    assert len(merged_single_track) == 1, (
        f"'SkipMultiTrackCells' should drop cell {CELL_A} and keep only cell {CELL_B}, i.e. 1 hit, "
        f"got {len(merged_single_track)}"
    )
    check_hit(OUT_COLL_SINGLE_TRACK, 0, merged_single_track[0], CELL_B, 4.0, 0.004, 1, 3.0)

    # --- RepresentativeKinematics = "Average" ------------------------------
    # Same sums as "SumAll", but the kinematics are now the unweighted mean over the summed hits, so
    # cell 1 sits at the mean of the times 2.0, 0.5, 1.0 and 1.5 ns, i.e. at 1.25 ns. Cell 2 has a
    # single hit and is therefore unchanged.
    merged_average = frame.get(OUT_COLL_AVERAGE)
    assert len(merged_average) == 2, (
        f"'Average' should give one hit per cell, i.e. 2, got {len(merged_average)}"
    )
    check_hit(OUT_COLL_AVERAGE, 0, merged_average[0], CELL_A, 3.75, 0.00375, 0, 1.25)
    check_hit(OUT_COLL_AVERAGE, 1, merged_average[1], CELL_B, 4.0, 0.004, 1, 3.0)

    print("[check] All assertions passed.")
    for coll_name in (
        OUT_COLL_ALL,
        OUT_COLL_PRIMARY,
        OUT_COLL_PER_TRACK,
        OUT_COLL_SINGLE_TRACK,
        OUT_COLL_AVERAGE,
    ):
        print(f"        {coll_name}: {len(frame.get(coll_name))} hits")


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
