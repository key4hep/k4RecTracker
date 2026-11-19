# Checks of the CombinatorialTrackFinder test output (test_CombinatorialTrackFinder.sh):
#   - all events were processed
#   - at least 8 of the 10 simulated muons have a track: >= 3 hits, >= 75% of them from the muon
#   - the median |pT(track) - pT(muon)| / pT(muon) of these tracks is below 5%
#   - every track has its AtFirstHit, AtLastHit and AtIP states
# usage: python3 checkCombinatorialTrackFinderTest.py <tracks.edm4hep.root> <compact.xml>
import math
import sys
import ROOT
from podio import root_io
import edm4hep

file_name = sys.argv[1]
compact_file = sys.argv[2]
expected_events = 10

# the field inside the tracker, from the geometry the events were simulated with
ROOT.gSystem.Load("libDDCore")
detector = ROOT.dd4hep.Detector.getInstance()
detector.fromXML(compact_file)
field_tesla = abs(
    detector.field().magneticField(ROOT.dd4hep.Position(0, 0, 0)).z() / ROOT.dd4hep.tesla
)
print(f"field at the IP: {field_tesla} T")

number_of_events = 0
muons_with_track = 0
relative_pt_errors = []
tracks_missing_states = 0
for event in root_io.Reader(file_name).get("events"):
    number_of_events += 1

    # the generated muon of this event
    muon = None
    for particle in event.get("MCParticles"):
        if particle.getGeneratorStatus() == 1 and abs(particle.getPDG()) == 13:
            muon = particle
    if muon is None:
        continue
    muon_momentum = muon.getMomentum()
    muon_pt = math.hypot(muon_momentum.x, muon_momentum.y)

    # reco hit (collection ID, index) -> MC particle that made it
    particle_of_hit = {}
    for link in event.get("InnerTrackerTrackHitSimLinks"):
        hit_id = link.getFrom().getObjectID()
        particle_of_hit[(hit_id.collectionID, hit_id.index)] = link.getTo().getParticle()

    found_track_of_muon = False
    for track in event.get("InnerTrackerTracks"):
        locations = [state.location for state in track.getTrackStates()]
        for wanted in (
            edm4hep.TrackState.AtFirstHit,
            edm4hep.TrackState.AtLastHit,
            edm4hep.TrackState.AtIP,
        ):
            if wanted not in locations:
                tracks_missing_states += 1
                break

        hits = list(track.getTrackerHits())
        hits_from_muon = 0
        for hit in hits:
            hit_id = hit.getObjectID()
            particle = particle_of_hit.get((hit_id.collectionID, hit_id.index))
            if particle is not None and particle.getObjectID() == muon.getObjectID():
                hits_from_muon += 1
        if len(hits) < 3 or hits_from_muon < 0.75 * len(hits) or found_track_of_muon:
            continue
        found_track_of_muon = True

        # pT from the curvature of the state at the first hit
        for state in track.getTrackStates():
            if state.location == edm4hep.TrackState.AtFirstHit and state.omega != 0:
                track_pt = 0.299792458 * field_tesla / (abs(state.omega) * 1000.0)
                relative_pt_errors.append(abs(track_pt - muon_pt) / muon_pt)
    if found_track_of_muon:
        muons_with_track += 1

relative_pt_errors.sort()
median_pt_error = relative_pt_errors[len(relative_pt_errors) // 2] if relative_pt_errors else 1.0
print(
    f"events: {number_of_events}, muons with a track: {muons_with_track}, "
    f"median |dpT|/pT: {median_pt_error:.4f}, tracks missing a state: {tracks_missing_states}"
)

failures = []
if number_of_events != expected_events:
    failures.append(f"expected {expected_events} events, found {number_of_events}")
if muons_with_track < 8:
    failures.append(f"only {muons_with_track} of {number_of_events} muons have a track")
if median_pt_error > 0.05:
    failures.append(f"median |dpT|/pT = {median_pt_error:.3f} is above 0.05")
if tracks_missing_states > 0:
    failures.append(f"{tracks_missing_states} tracks miss the AtFirstHit, AtLastHit or AtIP state")
for failure in failures:
    print("FAILED:", failure)
sys.exit(1 if failures else 0)
