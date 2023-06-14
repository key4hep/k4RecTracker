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

Cloning:

```bash
git clone https://github.com/mjbasso/k4RecTracker.git
```

Installing:

```bash
cd k4RecTracker
make
```

## Repository content

* `ARCdigi`: ARC digitization (for now, this step produces 'reco' collection)
* `VTXdigi`: vertex detector digitization (for now, this step produces 'reco' collection)
* `Tracking`: tracking algorithms orchestrating [GenFit](https://github.com/GenFit/GenFit)

## Execute Examples 

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
