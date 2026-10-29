# Environment for running the TrackingTest steering files.
# Source this AFTER setting up a key4hep environment (K4GEO must be defined).
# It puts the local k4RecTracker (GGTFTrackFinder, GenfitTrackFitter) and
# k4Reco (ConformalTracking, DDPlanarDigi, RefitFinal) builds first in the
# respective paths, so they take precedence over the stack copies.

K4RECTRACKER=/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker
K4RECO_LOCAL=/afs/cern.ch/work/a/adevita/public/workDir/k4Reco

export LD_LIBRARY_PATH=${K4RECTRACKER}/build/Tracking:${K4RECO_LOCAL}/build/k4Reco:${LD_LIBRARY_PATH}
export GAUDI_PLUGIN_PATH=${K4RECTRACKER}/build/Tracking:${K4RECO_LOCAL}/build/k4Reco:${GAUDI_PLUGIN_PATH}
export PYTHONPATH=${K4RECTRACKER}/build/Tracking/genConfDir:${K4RECTRACKER}/build/WireTrackerDigi/genConfDir:${K4RECO_LOCAL}/python:${K4RECO_LOCAL}/build/k4Reco/genConfDir:${PYTHONPATH}
