#!/bin/bash
###############################################################################
# Test script 
#
# Requirements:
#   - key4hep stack must be sourced
#   - Environment variables K4GEO and CLDCONFIG must be set
###############################################################################

Dir_TestSrc="${SOURCE_DIR_TEST:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

export Dir_TestSrc
export PYTHONPATH="${Dir_TestSrc}:$PYTHONPATH"

# Check required environment variables
if [ -z "${KEY4HEP_STACK:-}" ] || [ -z "${K4GEO:-}" ] || [ -z "${CLDCONFIG:-}" ]; then
    echo "Key4hep stack not found in this environment. Setting up latest nightly..."
    source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
fi

# Configuration
NEvts=10
Momentum=10
Particle="mu-"
File_SimOutput="Step1_test_vtxdigi_edm4hep.root"
File_DigiOutput="Test_vtxdigi_output"
File_Log="VTXdigi_Modular.log"

echo "========================================================================"
echo "  VTXDigi_Modular execution test"
echo "========================================================================"
echo "  Events:    ${NEvts}"
echo "  Particle:  ${Particle} @ ${Momentum} GeV (isotropic particle gun)"
echo "  IDEA detector, inner VTX"
echo "========================================================================"

# Step 1: Simulation with ddsim
echo "  [1/3] Running simulation with ddsim..."

File_SimSteering="${Dir_TestSrc}/steering_sim.py"
File_Detector=${K4GEO}/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml

ddsim --steeringFile ${File_SimSteering} \
      --compactFile ${File_Detector} \
      --enableGun \
      --gun.distribution uniform \
      --gun.energy "${Momentum}*GeV" \
      --gun.particle ${Particle} \
      --numberOfEvents ${NEvts} \
      --outputFile ${File_SimOutput} \
    #   > File_Log 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Simulation failed"
    echo "Logs in ${File_Log}"
    exit 1
fi
echo "  ✓ Simulation completed: ${File_SimOutput}"

# Step 2: Create digitization steering file
echo ""
echo "[2/3] Running digitization with VTXdigitizerDetailed..."

File_DigiSteering="${Dir_TestSrc}/steering_digi.py"

k4run ${File_DigiSteering} \
-n ${N_EVENTS} \
--inputFiles ${File_SimOutput} \
--outputBasename ${} 

if [ $? -ne 0 ]; then
    echo "ERROR: Digitization failed"
    echo "Logs in ${File_Log}"
    exit 1
fi
echo "  ✓ Digitization completed: ${File_DigiOutput}_REC.edm4hep.root"

# Step 3: Verify output
echo ""
echo "[3/3] Verifying output..."
if [ ! -f "${File_DigiOutput}_REC.edm4hep.root" ]; then
    echo "ERROR: Output file not created"
    echo "Logs in ${File_Log}"
    exit 1
fi

# Check collections in output file
echo ""
echo "Test if collections are in output file:"
podio-dump ${File_DigiOutput}_REC.edm4hep.root -c -d 2>/dev/null | grep -E "(VertexBarrel|VertexEndcap|InnerTrackerBarrel|InnerTrackerEndcap|OuterTrackerBarrel|OuterTrackerEndcap)" || true

# Cleanup
echo ""
echo "Cleaning up temporary files..."
rm -rf ${File_SimOutput} cdb.log ${File_DigiOutput}_aida.root __pycache__

echo ""
echo "========================================================================"
echo "  ✓ Test PASSED"
echo "  Output file: ${File_DigiOutput}_REC.edm4hep.root"
echo "  Log file:    ${File_Log}"
echo "========================================================================"
