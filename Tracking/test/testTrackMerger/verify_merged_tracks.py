import argparse
from podio import root_io

parser = argparse.ArgumentParser(description="Count tracks per event in an EDM4hep file.")
parser.add_argument("input_file", help="Input EDM4hep ROOT file")
parser.add_argument(
    "--collection",
    default="MarlinTrkTracks",
    help="Name of the track collection (default: MarlinTrkTracks)",
)

args = parser.parse_args()

reader = root_io.Reader(args.input_file)

for i, event in enumerate(reader.get("events")):
    tracks = event.get(args.collection)
    print(f"Event {i}: {len(tracks)} tracks")
    if len(tracks) == 0:
        raise ValueError("No entries found in track collection!")
