#!/bin/bash
# End-to-end test of CombinatorialTrackFinder: 10 simulated muons in CLD, digitised vertex,
# inner and outer tracker hits, track finding (runCombinatorialTrackFinder.py), then checks of the tracks.
set -e

# the steering files are next to this script; the outputs go to the current (build) directory
TEST_DIR=$(dirname "$0")
COMPACT=$K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml

# clean up previous output files
rm -f ctf_test_sim.edm4hep.root ctf_test_tracks.edm4hep.root

# one 10 GeV muon per event, 20 to 160 degrees (ddsim angles are in radians)
ddsim --compactFile $COMPACT --enableGun --gun.particle mu- --gun.energy "10*GeV" \
  --gun.distribution uniform --gun.thetaMin 0.349 --gun.thetaMax 2.793 \
  --numberOfEvents 10 --random.seed 42 --outputFile ctf_test_sim.edm4hep.root

k4run $TEST_DIR/runCombinatorialTrackFinder.py --inputFile ctf_test_sim.edm4hep.root \
  --outputFile ctf_test_tracks.edm4hep.root --digitise

python3 $TEST_DIR/checkCombinatorialTrackFinderTest.py ctf_test_tracks.edm4hep.root $COMPACT
