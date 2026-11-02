#!/bin/bash

# clean up previous output
rm -f testFitter.root

k4run runTestTrackFitter.py --inputFile ../testTrackFinder/out_tracks.root --outputFile testFitter.root
