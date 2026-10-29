#!/bin/bash
# Runs the Z -> uds tracking comparison for IDEA_o1_v04, following the
# production pipeline runSequence.sh (train branch), for NEV=1 event:
#   pythia (Zcard.cmd, e+e- @ 91 GeV) -> ddsim (crossing angle + beam-spot
#   smearing) -> digi -> {GGTF+Genfit, ConformalTracking+Refit} -> comparison
# Requires a key4hep environment to be set up beforehand (K4GEO must be defined).
set -euo pipefail
cd "$(dirname "$0")"

source setup_env.sh

NEV=1
SEED=1
GEN=gen_zuds
SIM=sim_zuds.edm4hep.root
DIGI=digi_zuds.edm4hep.root
AOUT=tracks_ggtf_genfit_zuds.edm4hep.root
BOUT=tracks_ct_refit_zuds.edm4hep.root

echo "=== 1/7 pythia: Z -> uds @ 91 GeV ($NEV event, seed $SEED) ==="
k4run pythia.py \
      -n $NEV \
      --IOSvc.Output ${GEN}.root \
      --HepMCFileWriter.Filename ${GEN}.hepmc \
      --Pythia8.PythiaInterface.pythiacard Zcard.cmd

echo "=== 2/7 ddsim ($NEV event) ==="
ddsim --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v04/IDEA_o1_v04.xml \
      --outputFile $SIM \
      --inputFiles ${GEN}.hepmc \
      --numberOfEvents $NEV \
      --vertexSigma 5.96e-3 23.8e-6 0.397 36.3 \
      --random.seed $SEED \
      --steeringFile SteeringFile_IDEA_o1_v04.py \
      --part.minimalKineticEnergy "0.00*MeV"
rm -f ${GEN}.root ${GEN}.hepmc

echo "=== 3/7 digitization (shared input) ==="
k4run runIDEA_digi.py --inputFile $SIM --outputFile $DIGI

echo "=== 4/7 pipeline A: GGTF + Genfit ==="
k4run runIDEA_GGTF_Genfit.py --inputFile $DIGI --outputFile $AOUT

echo "=== 5/7 pipeline B: ConformalTracking + Refit ==="
k4run runIDEA_ConformalTracking_Refit.py --inputFile $DIGI --outputFile $BOUT

echo "=== 6/7 comparison (stable charged truth) ==="
python3 compare_tracking.py --truth stable $AOUT $BOUT

echo "=== 7/7 plots ==="
python3 plot_tracking.py --truth stable --prefix tracking_zuds $AOUT $BOUT
