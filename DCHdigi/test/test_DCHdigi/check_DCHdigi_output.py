# file: check_DCHdigi_output.py
# author: Alvaro Tolosa-Delgado, CERN 2024
# to run: python3 check_DCHdigi_output.py
# goal: check distance hit-wire and print out a number:
#  0 : calculation of distance hit-wire was done properly
#  1 : problem with calculation of distance hit-wire
#  2 : problem with calculation of hit-projection position onto the wire

import ROOT

exit_code=0
# open debug output file generated by DCHdigi alg
f=ROOT.TFile("dch_digi_alg_debug.root")

# retrieve the hit-wire distance distribution
hDpw=f.Get("hDpw")
hDpw.Rebin(10)
# retrieve the distance at which the distance distribution has its maximum
distance_hit_wire_more_frequent=hDpw.GetXaxis().GetBinCenter( hDpw.GetMaximumBin() )
# the distance hit wire has the maximum around d=0.66cm
# if it is not the case, that means something weird is going on...
if 0.05 < abs(distance_hit_wire_more_frequent - 0.65) :
    exit_code+=1

# retrieve the hit-projection onto the wire to the wire distance distribution
hDww=f.Get("hDww")
# since the hit-projection onto the wire should be a point on the wire, the distance should be zero
# and all the counts are pushed to the bin number 1 (which is excluded from the integral)
hDww_integral=hDww.Integral(2,-1)
if hDww_integral != 0 :
    exit_code+=2

# we have to print the exit code, so it can be captured by the bash script
print(exit_code)