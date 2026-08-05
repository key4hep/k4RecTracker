#!/bin/bash

# clean up previous output files
rm -f k4RecTrackerTracking/test/testTrackFinder/out_sim_edm4hep.root
rm -f k4RecTrackerTracking/test/testTrackFinder/out_tracks.root

MODEL_PATH=$1

XML_FILE=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
STEERING_FILE=k4RecTrackerTracking/test/testTrackFinder/SteeringFile_IDEA_o1_v03.py
TBETA=0.6
TD=0.3

curl -o $STEERING_FILE https://raw.githubusercontent.com/key4hep/k4geo/master/example/SteeringFile_IDEA_o1_v03.py

ddsim --steeringFile $STEERING_FILE \
      --compactFile  $XML_FILE \
      -G --gun.distribution uniform --gun.particle mu- \
      --random.seed 42 \
      --numberOfEvents 1 \
      --outputFile k4RecTrackerTracking/test/testTrackFinder/out_sim_edm4hep.root

k4run k4RecTrackerTracking/test/testTrackFinder/runTestTrackFinder.py --inputFile k4RecTrackerTracking/test/testTrackFinder/out_sim_edm4hep.root --outputFile k4RecTrackerTracking/test/testTrackFinder/out_tracks.root --modelPath $MODEL_PATH --tbeta $TBETA --td $TD

