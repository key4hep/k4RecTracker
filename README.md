# k4RecTracker

This repository hosts Gaudi components related to vertex and tracker reconstruction as well as tracking.

## Dependencies

* ROOT
* PODIO
* EDM4HEP
* Gaudi
* k4FWCore
* DD4HEP

## Installation
```
#go somewhere public (to help support team helping you to debug your code)
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
git clone https://github.com/HEP-FCC/FCCDetectors
cd FCCDetectors
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ../../
export FCCDETECTORS=$PWD/FCCDetectors/;PATH=$PWD/FCCDetectors/install/bin/:$PATH;CMAKE_PREFIX_PATH=$PWD/FCCDetectors/install/:$CMAKE_PREFIX_PATH;LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib:$LD_LIBRARY_PATH;export PYTHONPATH=$PWD/FCCDetectors/install/python:$PYTHONPATH;LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib64:$LD_LIBRARY_PATH

git clone THIS REPO
cd k4RecTracker
# change branch to DChitGlobal
source setup.sh
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ../../
```

The following should be called in the folder hosting both `k4RecTracker` and `FCCDetectors` each time you start a new session (even you do not need to re-build): 
```
export FCCDETECTORS=$PWD/FCCDetectors/;PATH=$PWD/FCCDetectors/install/bin/:$PATH;CMAKE_PREFIX_PATH=$PWD/FCCDetectors/install/:$CMAKE_PREFIX_PATH;LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib:$LD_LIBRARY_PATH;export PYTHONPATH=$PWD/FCCDetectors/install/python:$PYTHONPATH;LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib64:$LD_LIBRARY_PATH
```

## Installation with makefile
(first copy locally FCCDetectors as above and export the relevant variables)

Cloning:

```bash
git clone git@github.com:key4hep/k4RecTracker.git
```

Installing:

```bash
cd k4RecTracker
make
```

Setting the environment:

```bash
source setup.sh
```

Fetching data:

```bash
make get_data
```

## Repository content

* `DCHdigi`: drift chamber digitization (for now, this step produces 'reco' collection)
* `ARCdigi`: ARC digitization (for now, this step produces 'reco' collection)
* `VTXdigi`: vertex detector digitization (for now, this step produces 'reco' collection)
* `Tracking`: tracking algorithms orchestrating [GenFit](https://github.com/GenFit/GenFit)

## Execute Examples 

```bash
k4run DCHdigi/test/runDCHsimpleDigitizer.py
```

```bash
k4run ARCdigi/test/runARCdigitizer.py
```

## Convention

For the syntax, try to follow the LLVM standards. You can format your code before to open a pull request with:

```bash
source /cvmfs/sft.cern.ch/lcg/contrib/clang/14.0.6/x86_64-centos7/setup.sh
clang-format -i path_to_your_file
```

## References:

These could perhaps be useful for newcomers:
1. [lhcb-98-064 COMP](https://cds.cern.ch/record/691746/files/lhcb-98-064.pdf)
2. [Hello World in the Gaudi Framework](https://lhcb.github.io/DevelopKit/02a-gaudi-helloworld)
