#!/bin/bash
set -uo pipefail
###############################################################################
# Test script
#
# Requirements:
#   - key4hep stack must be sourced
#   - Environment variable K4GEO must be set
#
# This test only checks if output collections exist, not whether they make any sense.
#
# TODO:
###############################################################################

Dir_TestSrc="${SOURCE_DIR_TEST:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

export Dir_TestSrc
export PYTHONPATH="${Dir_TestSrc}:${PYTHONPATH:-}"

# Check required environment variables
if [ -z "${KEY4HEP_STACK:-}" ] || [ -z "${K4GEO:-}" ]; then
    echo "Key4hep stack not found in this environment. Setting up latest nightly..."
    source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
fi

# Configuration
NEvts=10
Gun_Energy=10 # GeV
Gun_Particle="mu-"

File_SimOutput="Test_vtxdigi_M_simHits.root"
File_SimSteering="${Dir_TestSrc}/steering_sim.py"
File_Detector="${K4GEO}/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml"

File_DigiSteering="${Dir_TestSrc}/steering_digi.py"
File_LookupTable="${Dir_TestSrc}/lookup_table_dummy.init"
File_DigiOutput="Test_vtxdigi_M_output.root"

check_file() {
    local file_path="$1"
    if [ ! -f "$file_path" ]; then
        echo "ERROR: Required test input file not found: $file_path"
        exit 1
    fi
}

echo "========================================================================"
echo "  VTXDigi_Modular execution test"
echo "========================================================================"
echo "  Events:    ${NEvts}"
echo "  Particle:  ${Gun_Particle} @ ${Gun_Energy} GeV (particle gun uniform in cosTheta)"
echo "  ALLEGRO o1_v03, inner VTX"
echo "========================================================================"

echo "  Checking input files"

check_file "${File_SimSteering}"
check_file "${File_Detector}"
check_file "${File_DigiSteering}"
check_file "${File_LookupTable}"

rm -f "${File_SimOutput}"
rm -f "${File_DigiOutput}"


echo "  [1/3] Running simulation with ddsim..."

ddsim --steeringFile "${File_SimSteering}" \
    --compactFile "${File_Detector}" \
    --enableGun \
    --gun.distribution uniform \
    --gun.energy "${Gun_Energy}*GeV" \
    --gun.particle "${Gun_Particle}" \
    --numberOfEvents "${NEvts}" \
    --outputFile "${File_SimOutput}"

if [ $? -ne 0 ]; then
    echo "ERROR: Simulation failed."
    exit 1
fi

echo "  ✓ Simulation completed: ${File_SimOutput}"
echo ""



echo "  [2/3] Running digitization with VTXdigi_Modular using a lookup table..."

k4run "${File_DigiSteering}" \
    -n "${NEvts}" \
    --IOSvc.Input "${File_SimOutput}" \
    --IOSvc.Output "${File_DigiOutput}" \
    --GeoSvc.detectors "${File_Detector}" \
    --VTXBdigi_inner.LookupTableFile "${File_LookupTable}"

if [ $? -ne 0 ]; then
    echo "ERROR: Digitization failed."
    exit 1
fi

echo "  ✓ Digitization completed: ${File_DigiOutput}"



echo ""
echo "  [3/3] Verifying output..."
if [ ! -f "${File_DigiOutput}" ]; then
    echo "ERROR: Output file not created."
    exit 1
fi

## this check only verifies that the collections exists, but still passes on empty collections
# dump=$(podio-dump "${File_DigiOutput}") || {
#     echo "ERROR: podio-dump failed"
#     exit 1
# }
# echo "$dump" | grep -E "VTXBDigis|VTXBSimDigiLinks" || {
#     echo "ERROR: digitizer output collections missing."
#     exit 1
# }

## this checks for non-0 entry counts in the collections
python3 - "${File_DigiOutput}" <<'EOF' || exit 1
import sys
from podio.reading import get_reader

reader = get_reader(sys.argv[1])
n_events = n_sim = n_dig = n_lnk = 0
for frame in reader.get("events"):
    n_events += 1
    n_sim += len(frame.get("VertexBarrelCollection"))
    n_dig += len(frame.get("VTXBDigis"))
    n_lnk += len(frame.get("VTXBSimDigiLinks"))

print(f"  ✓ Found {n_events} events, {n_sim} simHits, {n_dig} digiHits, {n_lnk} links")

if n_events == 0:
    sys.exit("ERROR: no events in output file")
if n_sim == 0:
    sys.exit("ERROR: no sim hits in VertexBarrelCollection - simulation/geometry problem, not the digitizer")
if n_dig == 0 or n_lnk == 0:
    sys.exit("ERROR: sim hits present but digitizer produced no output")
EOF



echo ""
echo "  Cleaning up temporary files..."
rm -rf "${File_SimOutput}" __pycache__

echo ""
echo "========================================================================"
echo "  ✓ Test PASSED"
echo "  Output file: ${File_DigiOutput}"
echo "========================================================================"
