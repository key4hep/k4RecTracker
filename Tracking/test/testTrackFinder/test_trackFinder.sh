#!/bin/bash

XML_FILE=$K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
STEERING_FILE=IDEAsteeringFile.py
MODEL_PATH=/eos/project/f/fccsw-web/www/filesForSimDigiReco/IDEA/IDEA_o1_v03/SimpleGatrIDEAv3o1.onnx
TBETA=0.6
TD=0.3

k4run pythia.py -n 1 --Dumper.Filename out_pythia.hepmc --Pythia8.PythiaInterface.pythiacard Zcard.cmd

ddsim   --compactFile $XML_FILE \
        --outputFile out_sim_edm4hep.root \
        --inputFiles out_pythia.hepmc \
        --numberOfEvents 1 \
        --random.seed 42 \
        --part.minimalKineticEnergy "0.001*MeV" \
        --steeringFile  $STEERING_FILE
    
k4run runTestTrackFinder.py --inputFile out_sim_edm4hep.root --outputFile out_tracks.root --modelPath $MODEL_PATH --tbeta $TBETA --td $TD
