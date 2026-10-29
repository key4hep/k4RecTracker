"""
Pythia8, integrated in the FCCSW framework.

Generates according to a pythia .cmd file and saves them in fcc edm format.

"""

import os
from GaudiKernel import SystemOfUnits as units
from Gaudi.Configuration import *
from edm4hep import labels as e4_labels

from Configurables import EventDataSvc
from k4FWCore import ApplicationMgr, IOSvc

ApplicationMgr().EvtSel = "NONE"
ApplicationMgr().EvtMax = 2
ApplicationMgr().OutputLevel = INFO
ApplicationMgr().ExtSvc += ["RndmGenSvc", EventDataSvc("EventDataSvc")]

from Configurables import EventHeaderCreator

eventHeaderCreator = EventHeaderCreator("eventHeaderCreator", runNumber=42, eventNumberOffset=42)
ApplicationMgr().TopAlg += [eventHeaderCreator]


# from Configurables import GaussSmearVertex
# smeartool = GaussSmearVertex()
# smeartool.xVertexSigma = 5.96e-3 * units.mm
# smeartool.yVertexSigma = 23.8e-6 * units.mm
# smeartool.zVertexSigma = 0.397 * units.mm
# smeartool.tVertexSigma = 36.3 * units.picosecond

from Configurables import PythiaInterface

pythia8gentool = PythiaInterface()
### Example of pythia configuration file to generate events
# take from $K4GEN if defined, locally if not
path_to_pythiafile = os.environ.get("K4GEN", "")
pythiafilename = "Pythia_standard.cmd"
pythiafile = os.path.join(path_to_pythiafile, pythiafilename)
# Example of pythia configuration file to read LH event file
# pythiafile="options/Pythia_LHEinput.cmd"
pythia8gentool.pythiacard = pythiafile
pythia8gentool.doEvtGenDecays = False
pythia8gentool.printPythiaStatistics = True
pythia8gentool.pythiaExtraSettings = [""]

from Configurables import GenAlg

pythia8gen = GenAlg("Pythia8")
pythia8gen.SignalProvider = pythia8gentool
# pythia8gen.VertexSmearingTool = smeartool
pythia8gen.hepmc.Path = "hepmc"
ApplicationMgr().TopAlg += [pythia8gen]

# Store the generated event as HepMC3 as well as EDM4hep.  HepMC3 is used as
# the ddsim input because it preserves displaced production and decay vertices
# with Key4hep/DD4hep versions whose EDM4hep input reader collapses them.
from Configurables import HepMCFileWriter

hepmc_writer = HepMCFileWriter("HepMCFileWriter")
hepmc_writer.hepmc.Path = "hepmc"
hepmc_writer.Filename = "output_pythia.hepmc"
ApplicationMgr().TopAlg += [hepmc_writer]

### Reads an HepMC::GenEvent from the data service and writes a collection of EDM Particles
from Configurables import HepMCToEDMConverter

hepmc_converter = HepMCToEDMConverter()
hepmc_converter.hepmc.Path = "hepmc"
hepmc_converter.hepmcStatusList = []  # convert particles with all statuses
hepmc_converter.GenParticles.Path = e4_labels.MCParticles
ApplicationMgr().TopAlg += [hepmc_converter]

### Filters generated particles
# accept is a list of particle statuses that should be accepted
from Configurables import GenParticleFilter

genfilter = GenParticleFilter("StableParticles")
genfilter.accept = [1]
genfilter.GenParticles.Path = e4_labels.MCParticles
genfilter.GenParticlesFiltered.Path = "MCParticlesStable"
ApplicationMgr().TopAlg += [genfilter]

iosvc = IOSvc()
iosvc.Output = "output_pythia.root"
iosvc.outputCommands = ["keep *"]
