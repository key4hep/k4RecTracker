#!/bin/bash

# clean up previous output files
rm -rf CLDConfig  k4GaudiPandora cdb.log  output_pandora_bremRecovered.root

git clone https://github.com/key4hep/CLDConfig.git

git clone https://github.com/key4hep/k4GaudiPandora.git

cd CLDConfig/CLDConfig

ddsim --compactFile $K4GEO/FCCee/CLD/compact/CLD_o2_v07/CLD_o2_v07.xml --steeringFile cld_steer.py --enableGun --gun.particle e- --gun.energy 10*GeV --gun.direction 0.966,0.259,0.01 --numberOfEvents 10

k4run CLDReconstruction.py --inputFiles CLD_SIM.edm4hep.root --outputBasename CLD_SIM.edm4hep.Output --num-events -1

k4run ../../k4GaudiPandora/k4GaudiPandora/options/runPandora.py --IOSvc.Input=CLD_SIM.edm4hep.Output_REC.edm4hep.root --IOSvc.Output=output_pandora.root

cd ../../

k4run ElectronDecorator.py
