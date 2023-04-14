# k4RecTracker


This repository hosts Gaudi components related to vertex and tracker reconstruction as well as tracking.


## Dependencies

* ROOT
* PODIO
* EDM4HEP
* Gaudi
* k4FWCore

## Installation

```
source /cvmfs/sw.hsf.org/key4hep/setup.sh
git clone git@github.com:BrieucF/k4RecTracker.git
cd k4RecTracker
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 4
cd ..
LD_LIBRARY_PATH=$PWD/install/lib:$PWD/install/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/install/python:$PYTHONPATH

```
## Repository content

* `DCHdigi`: drift chamber digitization (for now, this step produces 'reco' collection)
* `VTXdigi`: vertex detector digitization (for now, this step produces 'reco' collection)
* `Tracking`: tracking algorithms orchestrating [GenFit](https://github.com/GenFit/GenFit)


## Execute Examples 


```
k4run DCHdigi/test/runDCHdigitizer.py

```

## References:
These could perhaps be usefule for newcomers. 
1. [lhcb-98-064 COMP](https://cds.cern.ch/record/691746/files/lhcb-98-064.pdf)
2. [Hello World in the Gaudi Framework](https://lhcb.github.io/DevelopKit/02a-gaudi-helloworld)
