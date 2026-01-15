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

set -e  # Exit on error
set +u  # Do not Exit on undefined variable while sourcing setup

# Save the initial directory
INITIAL_DIR=$(pwd)
# Get the directory of the script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Go to the script directory (to have a clean environment with access to local files)
cd "${SCRIPT_DIR}"

# Check required environment variables
if [ -z "${K4GEO:-}" ] || [ -z "${CLDCONFIG:-}" ]; then
    source /cvmfs/sw.hsf.org/key4hep/setup.sh
fi

set -u  # Exit on undefined variable

# Add script directory to PYTHONPATH to find local python modules
export PYTHONPATH="${SCRIPT_DIR}:$PYTHONPATH"

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

cat > collections_rec_level.txt << EOF
    BuildUpVertices                                    Vertex
     BuildUpVertices_RP                                 ReconstructedParticle
     BuildUpVertices_V0                                 Vertex
     BuildUpVertices_V0_RP                              ReconstructedParticle
     CalohitMCTruthLink                                 LCRelation[CalorimeterHit,MCParticle]
     ClusterMCTruthLink                                 LCRelation[Cluster,MCParticle]
     DebugHits                                          TrackerHitPlane
     ECalBarrelCollection                               SimCalorimeterHit
     ECALBarrel                                         CalorimeterHit
     ECALEndcap                                         CalorimeterHit
     ECalEndcapCollection                               SimCalorimeterHit
     EfficientMCParticles                               MCParticle
     HCalBarrelCollection                               SimCalorimeterHit
     HCALEndcap                                         CalorimeterHit
     HCALBarrel                                         CalorimeterHit
     HCALOther                                          CalorimeterHit
     HCalEndcapCollection                               SimCalorimeterHit
     HCalRingCollection                                 SimCalorimeterHit
     ITrackerEndcapHits                                 TrackerHitPlane
     ITrackerHits                                       TrackerHitPlane
     InefficientMCParticles                             MCParticle
     InnerTrackerBarrelCollection                       SimTrackerHit
     InnerTrackerBarrelHitsRelations                    LCRelation[TrackerHitPlane,SimTrackerHit]
     InnerTrackerEndcapCollection                       SimTrackerHit
     InnerTrackerEndcapHitsRelations                    LCRelation[TrackerHitPlane,SimTrackerHit]
     LooseSelectedPandoraPFOs                           ReconstructedParticle
     LumiCalClusters                                    Cluster
     LumiCalCollection                                  SimCalorimeterHit
     LumiCalHits                                        CalorimeterHit
     LumiCalRecoParticles                               ReconstructedParticle
     MCParticle                                         MCParticle
     MCParticlesSkimmed                                 MCParticle
     MCPhysicsParticles                                 MCParticle
     MCTruthClusterLink                                 LCRelation[MCParticle,Cluster]
     MCTruthRecoLink                                    LCRelation[MCParticle,ReconstructedParticle]
     MCTruthSiTracksLink                                LCRelation[MCParticle,Track]
     MUON                                               CalorimeterHit
     OTrackerEndcapHits                                 TrackerHitPlane
     OTrackerHits                                       TrackerHitPlane
     OuterTrackerBarrelCollection                       SimTrackerHit
     OuterTrackerBarrelHitsRelations                    LCRelation[TrackerHitPlane,SimTrackerHit]
     OuterTrackerEndcapCollection                       SimTrackerHit
     OuterTrackerEndcapHitsRelations                    LCRelation[TrackerHitPlane,SimTrackerHit]
     PandoraClusters                                    Cluster
     PandoraPFOs                                        ReconstructedParticle
     PandoraStartVertices                               Vertex
     PFOsFromJets                                       ReconstructedParticle
     PrimaryVertices                                    Vertex
     PrimaryVertices_RP                                 ReconstructedParticle
     RecoMCTruthLink                                    LCRelation[ReconstructedParticle,MCParticle]
     RefinedVertex                                      RefinedVertexJets|SingleVertexProbability
     RefinedVertexJets                                  ReconstructedParticle
     RefinedVertexJets_rel                              LCRelation[ReconstructedParticle,Vertex]
     RefinedVertexJets_vtx                              Vertex
     RefinedVertexJets_vtx_RP                           ReconstructedParticle
     RefinedVertices                                    Vertex
     RefinedVertices_RP                                 ReconstructedParticle
     RelationCaloHit                                    LCRelation[CalorimeterHit,SimCalorimeterHit]
     RelationMuonHit                                    LCRelation[CalorimeterHit,SimCalorimeterHit]
     SelectedPandoraPFOs                                ReconstructedParticle
     SiTracks                                           Track
     SiTracksCT                                         Track
     SiTracksMCTruthLink                                LCRelation[Track,MCParticle]
     SiTracks_Refitted                                  Track
     TightSelectedPandoraPFOs                           ReconstructedParticle
     TrackerHitPlane                                    TrackerHitPlane
     VXDTrackerHitRelations                             LCRelation[TrackerHitPlane,SimTrackerHit]
     VXDTrackerHits                                     TrackerHitPlane
     VXDEndcapTrackerHitRelations                       LCRelation[TrackerHitPlane,SimTrackerHit]
     VXDEndcapTrackerHits                               TrackerHitPlane
     VertexBarrelCollection                             SimTrackerHit
     VertexEndcapCollection                             SimTrackerHit
     VertexJets                                         ReconstructedParticle
     yth                                                RefinedVertexJets|y01,y12,y23,y34,y45,y56,y67,y78,y89,y910
     yth                                                VertexJets|y01,y12,y23,y34,y45,y56,y67,y78,y89,y910
     YokeBarrelCollection                               SimCalorimeterHit
     YokeEndcapCollection                               SimCalorimeterHit

EOF

cat > test_digi_steering.py << EOF
import os
from Gaudi.Configuration import INFO, WARNING, DEBUG

from Gaudi.Configurables import EventDataSvc, MarlinProcessorWrapper, GeoSvc, TrackingCellIDEncodingSvc
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from py_utils import SequenceLoader, parse_collection_patch_file
from k4MarlinWrapper.io_helpers import IOHandlerHelper

parser_group = parser.add_argument_group("CLDReconstruction.py custom options")
# Need the dummy input such that the IOHandlerHelper.add_reader call below does not crash when called with --help
parser_group.add_argument("--inputFiles", action="store", nargs="+", metavar=("file1", "file2"), help="One or multiple input files", default=["dummy_input.edm4hep.root"])
parser_group.add_argument("--outputBasename", help="Basename of the output file(s)", default="output")
parser_group.add_argument("--compactFile", help="Compact detector file to use", type=str, default=os.environ["K4GEO"] + "/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml")
tracking_group = parser_group.add_mutually_exclusive_group()
reco_args = parser.parse_known_args()[0]

evtsvc = EventDataSvc("EventDataSvc")
iosvc = IOSvc()

svcList = [evtsvc, iosvc]
algList = []

CONFIG = {
             "CalorimeterIntegrationTimeWindow": "10ns",
             "CalorimeterIntegrationTimeWindowChoices": ["10ns", "400ns"],
             "Overlay": "False",
             "OverlayChoices": ["False", "91GeV", "365GeV"],
             "VertexUnconstrained": "OFF",
             "VertexUnconstrainedChoices": ["ON", "OFF"],
             "OutputMode": "EDM4Hep",
             "OutputModeChoices": ["LCIO", "EDM4hep"] #, "both"] FIXME: both is not implemented yet
}

REC_COLLECTION_CONTENTS_FILE = "collections_rec_level.txt" # file with the collections to be patched in when writing from LCIO to EDM4hep

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [reco_args.compactFile]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False
svcList.append(geoservice)

cellIDSvc = TrackingCellIDEncodingSvc("CellIDSvc")
cellIDSvc.EncodingStringParameterName = "GlobalTrackerReadoutID"
cellIDSvc.GeoSvcName = geoservice.name()
cellIDSvc.OutputLevel = INFO
svcList.append(cellIDSvc)

if len(geoservice.detectors) > 1:
    # we are making assumptions for reconstruction parameters based on the detector option, so we limit the possibilities
    raise RuntimeError("Too many XML files for the detector path, please only specify the main file!")

# from https://github.com/HEP-FCC/FCCeePhysicsPerformance/blob/d6ecee2c2c3ed5d76db55a3ae18ced349b2b914a/General/README.md?plain=1#L457-L467
# for december 2022

sequenceLoader = SequenceLoader(
    algList,
    # global_vars can be used in sequence-loaded modules without explicit import
    global_vars={"CONFIG": CONFIG, "geoservice": geoservice, "reco_args": reco_args,},
)

io_handler = IOHandlerHelper(algList, iosvc)
io_handler.add_reader(reco_args.inputFiles)

MyAIDAProcessor = MarlinProcessorWrapper("MyAIDAProcessor")
MyAIDAProcessor.OutputLevel = WARNING
MyAIDAProcessor.ProcessorType = "AIDAProcessor"
MyAIDAProcessor.Parameters = {
                              "Compress": ["1"],
                              "FileName": [f"{reco_args.outputBasename}_aida"],
                              "FileType": ["root"]
                              }

EventNumber = MarlinProcessorWrapper("EventNumber")
EventNumber.OutputLevel = WARNING
EventNumber.ProcessorType = "Statusmonitor"
EventNumber.Parameters = {
                          "HowOften": ["1"]
                          }

# setup AIDA histogramming and add eventual background overlay
algList.append(MyAIDAProcessor)
sequenceLoader.load("Overlay")
# tracker hit digitisation
sequenceLoader.load("TrackingDigi")

# event number processor, down here to attach the conversion back to edm4hep to it
algList.append(EventNumber)

DST_KEEPLIST = ["MCParticlesSkimmed", "MCPhysicsParticles"]

DST_SUBSETLIST = ["EfficientMCParticles", "InefficientMCParticles", "MCPhysicsParticles"]

# Make sure that all collections are always available by patching in missing ones on-the-fly
collPatcherRec = MarlinProcessorWrapper(
    "CollPatcherREC", OutputLevel=INFO, ProcessorType="PatchCollections"
)
collPatcherRec.Parameters = {
    "PatchCollections": parse_collection_patch_file(REC_COLLECTION_CONTENTS_FILE)
}
algList.append(collPatcherRec)

io_handler.add_edm4hep_writer(f"{reco_args.outputBasename}_REC.edm4hep.root", ["keep *"])
# <DST output for edm4hep>


# We need to attach all the necessary converters
io_handler.finalize_converters()

ApplicationMgr( TopAlg = algList,
                EvtSel = 'NONE',
                EvtMax = 3, # Overridden by the --num-events switch to k4run
                ExtSvc = svcList,
                OutputLevel=WARNING
              )
EOF

# Run digitization
k4run test_digi_steering.py --inputFiles ${TEMP_SIM} --outputBasename ${OUTPUT_FILE} -n ${N_EVENTS}

if [ $? -ne 0 ]; then
    echo "ERROR: Digitization failed"
    rm -f ${TEMP_SIM} test_digi_steering.py
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
rm -rf ${TEMP_SIM} test_digi_steering.py cdb.log ${OUTPUT_FILE}_aida.root collections_rec_level.txt __pycache__

echo ""
echo "========================================================================"
echo "  ✓ Test PASSED"
echo "  Output file: ${OUTPUT_FILE}_REC.edm4hep.root"
echo "========================================================================"

# If use source to lauch the script, avoid exiting the parent shell on error or undefined variable
set +e
set +u

# Return to the initial directory
cd "${INITIAL_DIR}"