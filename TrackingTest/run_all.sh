#!/bin/bash
# Runs the muon-gun tracking comparison for IDEA_o1_v04:
#   sim (10 events x 3 mu- @ 10 GeV) -> digi -> {GGTF+Genfit, ConformalTracking+Refit} -> comparison
# Requires a key4hep environment to be set up beforehand (K4GEO must be defined).
set -euo pipefail
cd "$(dirname "$0")"

source setup_env.sh

SIM=sim.edm4hep.root
DIGI=digi.edm4hep.root
AOUT=tracks_ggtf_genfit.edm4hep.root
BOUT=tracks_ct_refit.edm4hep.root

echo "=== 1/6 ddsim: 10 events, 3 mu- @ 10 GeV ==="
ddsim --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v04/IDEA_o1_v04.xml \
      -G --numberOfEvents 10 --outputFile $SIM \
      --gun.energy '10*GeV' --gun.particle 'mu-' --gun.multiplicity 3 \
      --gun.distribution uniform --random.seed 42 \
      --steeringFile ddsimSteeringIDEA.py

echo "=== 2/6 digitization (shared input) ==="
k4run runIDEA_digi.py --inputFile $SIM --outputFile $DIGI

echo "=== 3/6 pipeline A: GGTF + Genfit ==="
k4run runIDEA_GGTF_Genfit.py --inputFile $DIGI --outputFile $AOUT

echo "=== 4/6 pipeline B: ConformalTracking + Refit ==="
k4run runIDEA_ConformalTracking_Refit.py --inputFile $DIGI --outputFile $BOUT

echo "=== 5/6 comparison ==="
python3 compare_tracking.py $AOUT $BOUT

echo "=== 6/6 plots ==="
python3 plot_tracking.py $AOUT $BOUT
