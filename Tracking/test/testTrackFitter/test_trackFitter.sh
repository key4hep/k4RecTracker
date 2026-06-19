#!/bin/bash

# clean up old temp files
rm -f testFitter.root

k4run runTestTrackFitter.py --inputFile ../testTrackFinder/out_tracks.root