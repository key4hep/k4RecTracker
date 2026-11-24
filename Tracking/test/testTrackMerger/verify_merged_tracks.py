import argparse
from collections import Counter

import dd4hep
from podio.reading import get_reader

CELLID_ENCODING = "system:5,side:-2,layer:9,module:8,sensor:8"

TRACK_COLLECTION = "MarlinTrkTracks"

# These thresholds are tuned to be relatively loose and only fire if things
# break entirely

# Minimum number of tracks with silicon and TPC hits in every event
MIN_SI_TPC_TRACKS_PER_EVENT = 2
# Minimum average number of tracks per event
MIN_MEAN_TRACKS = 8.0
# Minimum fraction of all tracks (over all events) with silicon and TPC hits
MIN_SI_TPC_FRACTION = 0.80

# Sub-detector system IDs as defined in k4geo FCCee/ILD_FCCee/compact/ILD_common_FCCee/basic_defs.xml.
SYSTEM_NAMES = {
    1: "VXDBarrel",
    2: "VXDEndcap",
    3: "ITBarrel",
    4: "TPC",
    5: "ITEndcap",
    6: "SET",
}
SILICON_SYSTEMS = {1, 2, 3, 5}
TPC_SYSTEMS = {4}


def system_name(system):
    return SYSTEM_NAMES.get(system, f"unknown({system})")


def hit_systems(track, decoder):
    """Get the set of sub-detector system IDs of all hits of a track"""
    return {decoder.get(hit.getCellID(), "system") for hit in track.getTrackerHits()}


def is_si_tpc_track(systems):
    """Whether a track has hits in the silicon trackers as well as in the TPC"""
    return bool(systems & SILICON_SYSTEMS) and bool(systems & TPC_SYSTEMS)


parser = argparse.ArgumentParser(
    description="Check that the tracks in an EDM4hep file contain merged silicon + TPC tracks."
)
parser.add_argument("input_file", help="Input EDM4hep ROOT file")
parser.add_argument(
    "--verbose",
    "-v",
    action="store_true",
    help="Print the sub-detector combinations of the tracks for every event, not only the total",
)
args = parser.parse_args()

decoder = dd4hep.core.DDSegmentation.BitFieldCoder(CELLID_ENCODING)

reader = get_reader(args.input_file)

errors = []
n_tracks_total = 0
n_si_tpc_total = 0
combinations_total = Counter()

n_events = 0
events = reader.get("events")
for i, event in enumerate(events):
    tracks = event.get(TRACK_COLLECTION)
    n_si_tpc = 0
    combinations = Counter()

    if len(tracks) == 0:
        errors.append(f"Event {i}: no tracks in collection '{TRACK_COLLECTION}'")

    for track in tracks:
        systems = hit_systems(track, decoder)
        combinations["+".join(system_name(s) for s in sorted(systems)) or "no hits"] += 1
        if is_si_tpc_track(systems):
            n_si_tpc += 1

    if args.verbose:
        print(f"Event {i}: {len(tracks)} tracks, {n_si_tpc} with silicon and TPC hits")
        for combination, count in combinations.most_common():
            print(f"    {count:3d} x {combination}")

    if n_si_tpc < MIN_SI_TPC_TRACKS_PER_EVENT:
        errors.append(
            f"Event {i}: only {n_si_tpc} tracks with silicon and TPC hits "
            f"(expected at least {MIN_SI_TPC_TRACKS_PER_EVENT})"
        )

    n_tracks_total += len(tracks)
    n_si_tpc_total += n_si_tpc
    combinations_total.update(combinations)

n_events = len(events)
mean_tracks = n_tracks_total / n_events if n_events > 0 else 0.0
si_tpc_fraction = n_si_tpc_total / n_tracks_total if n_tracks_total > 0 else 0.0
print(
    f"Total: {n_events} events, {n_tracks_total} tracks ({mean_tracks:.2f} / event), "
    f"{n_si_tpc_total} ({si_tpc_fraction:.1%}) with silicon and TPC hits"
)
for combination, count in combinations_total.most_common():
    print(f"    {count:3d} x {combination}")

if n_events == 0:
    errors.append("No events found in input file")

if mean_tracks < MIN_MEAN_TRACKS:
    errors.append(
        f"Only {mean_tracks:.2f} tracks per event on average (expected at least {MIN_MEAN_TRACKS:.2f})"
    )

if si_tpc_fraction < MIN_SI_TPC_FRACTION:
    errors.append(
        f"Only {si_tpc_fraction:.1%} of all tracks have silicon and TPC hits "
        f"(expected at least {MIN_SI_TPC_FRACTION:.1%})"
    )

if errors:
    raise ValueError("Verification of merged tracks failed:\n" + "\n".join(errors))
