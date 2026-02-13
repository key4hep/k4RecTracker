#!/bin/bash
###############################################################################
# Test script for VTXdigitizerDetailed
#
# This script runs a quick test of the VTXdigitizerDetailed digitizer
# Usage: ./test_VTXdigitizerDetailed.sh
#
# Requirements:
#   - key4hep stack must be sourced
#   - Environment variables K4GEO and CLDCONFIG must be set
###############################################################################


SRC_TEST="${SOURCE_DIR_TEST:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export SRC_TEST
export PYTHONPATH="${SRC_TEST}:$PYTHONPATH"

# Check required environment variables
if [ -z "${K4GEO:-}" ] || [ -z "${CLDCONFIG:-}" ]; then
    source /cvmfs/sw.hsf.org/key4hep/setup.sh
fi

# Configuration
N_EVENTS=10
MOMENTUM=10
PARTICLE="mu-"
TEMP_SIM="Step1_test_vtxdigi_edm4hep.root"
OUTPUT_FILE="Test_vtxdigi_output"

echo "========================================================================"
echo "  VTXdigitizerDetailed Quick Test"
echo "========================================================================"
echo "  Events:    ${N_EVENTS}"
echo "  Particle:  ${PARTICLE} @ ${MOMENTUM} GeV"
echo "========================================================================"

# Step 1: Simulation with ddsim
echo ""
echo "[1/3] Running simulation with ddsim..."
ddsim --steeringFile ${CLDCONFIG}/share/CLDConfig/cld_steer.py \
      --compactFile ${K4GEO}/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml \
      --enableGun \
      --gun.distribution uniform \
      --gun.energy "${MOMENTUM}*GeV" \
      --gun.particle ${PARTICLE} \
      --numberOfEvents ${N_EVENTS} \
      --outputFile ${TEMP_SIM} \
      > /dev/null 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Simulation failed"
    exit 1
fi
echo "  ✓ Simulation completed: ${TEMP_SIM}"

# Step 2: Create digitization steering file
echo ""
echo "[2/3] Running digitization with VTXdigitizerDetailed..."

# Run digitization
k4run ${SRC_TEST}/test_digi_steering.py --inputFiles ${TEMP_SIM} --outputBasename ${OUTPUT_FILE} -n ${N_EVENTS}

if [ $? -ne 0 ]; then
    echo "ERROR: Digitization failed"
    exit 1
fi
echo "  ✓ Digitization completed: ${OUTPUT_FILE}_REC.edm4hep.root"

# Step 3: Verify output
echo ""
echo "[3/3] Verifying output..."
if [ ! -f "${OUTPUT_FILE}_REC.edm4hep.root" ]; then
    echo "ERROR: Output file not created"
    exit 1
fi

# Check collections in output file
echo ""
echo "Test if collections are in output file:"
podio-dump ${OUTPUT_FILE}_REC.edm4hep.root -c -d 2>/dev/null | grep -E "(VertexBarrel|VertexEndcap|InnerTrackerBarrel|InnerTrackerEndcap|OuterTrackerBarrel|OuterTrackerEndcap)" || true

# Cleanup
echo ""
echo "Cleaning up temporary files..."
rm -rf ${TEMP_SIM} cdb.log ${OUTPUT_FILE}_aida.root __pycache__

echo ""
echo "========================================================================"
echo "  ✓ Test PASSED"
echo "  Output file: ${OUTPUT_FILE}_REC.edm4hep.root"
echo "========================================================================"
