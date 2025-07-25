#!/bin/bash

XML_FILE=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
STEERING_FILE=SteeringFile_IDEA_o1_v03.py
MODEL_PATH=/eos/project/f/fccsw-web/www/filesForSimDigiReco/IDEA/IDEA_o1_v03/SimpleGatrIDEAv3o1.onnx
TBETA=0.6
TD=0.3

ddsim --steeringFile $STEERING_FILE \
      --compactFile  $XML_FILE \
      -G --gun.distribution uniform --gun.particle e- \
      --random.seed 42 \
      --numberOfEvents 1 \
      --outputFile out_sim_edm4hep.root 
    
k4run runTestTrackFinder.py --inputFile out_sim_edm4hep.root --outputFile out_tracks.root --modelPath $MODEL_PATH --tbeta $TBETA --td $TD
