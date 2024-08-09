# file: test_DCHdigi.sh
# author: Alvaro Tolosa-Delgado, CERN 2024
# to run: sh + test_DCHdigi.sh
# goal: run sim-digitizer of the DCH v2, and return code printed by check_DCHdigi_output.py

# run simulation with the drift chamber alone
ddsim --steeringFile sim_steering.py --outputFile 'dch_proton_10GeV.root' -N 10 --runType batch --random.seed 42

# download file for cluster counting technique
wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root

# run digitizer for position smearing and cluster counting calculation
k4run runDCHdigi.py

# check distribution of distance from hit position to the wire
check_DCHdigi_output=$( (python3 check_DCHdigi_output.py) 2>&1)

# return value printed out by the previous python script
exit $check_DCHdigi_output
