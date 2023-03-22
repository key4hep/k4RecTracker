# k4-project-template


This repository can be a starting point and template for projects using the Key4HEP software stack.



## Dependencies

* ROOT

* PODIO

* Gaudi

* k4FWCore

## Installation


```
source /cvmfs/sw.hsf.org/key4hep/setup.sh
mkdir build install
cd build;
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install


```

## Execute Examples 


```
k4run ../k4TestFWCore/options/createHelloWorld.py 

k4run ../k4TestFWCore/options/createExampleEventData.py 

```


## References:
These could perhaps be usefule for newcomers. 
1. [lhcb-98-064 COMP](https://cds.cern.ch/record/691746/files/lhcb-98-064.pdf)
2. [Hello World in the Gaudi Framework](https://lhcb.github.io/DevelopKit/02a-gaudi-helloworld)
