import sys
import math
import ROOT
from array import array
from ROOT import TFile, TTree
import numpy as np
from podio import root_io
import edm4hep

print("started reading")
input_file = "/eos/experiment/fcc/ee/datasets/DC_tracking/Pythia_evaluation/scratch/Zcard_fakeCalo_test/1/output_tracking.root"
reader = root_io.Reader(input_file)
print("finished reading")
coordinates_tracks = []
for i, event in enumerate(reader.get("events")):
    print("event", i)
    tracks = event.get("CDCHTracks")
    for j, track in enumerate(tracks):
        print("track", j)
        hits_in_track = track.getTrackerHits()
        coordinates = []
        for hit in hits_in_track:
            position = hit.getPosition()
            x = position.x
            y = position.y
            z = position.z
            coordinates.append([x, y, z])
        coordinates = np.array(coordinates)
        np.save("track" + str(j) + ".npy", coordinates)
