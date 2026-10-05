from DDSim.DD4hepSimulation import DD4hepSimulation

###################################
# user options
simulateThinSiliconDetectors = ["VertexBarrel"]
simulateThinSilicon = True  # False -> don't use thin silicon settings in barrel
###################################

SIM = DD4hepSimulation()
SIM.printLevel = 3
SIM.runType = "batch"

## Lorentz boost for the crossing angle, in radian!
SIM.crossingAngleBoost = 0

SIM.enableDetailedShowerMode = False
SIM.enableG4GPS = False
SIM.enableG4Gun = False
SIM.enableGun = False

## Physics list to use in simulation


if not simulateThinSilicon:
    SIM.physicsList = "FTFP_BERT"
else:
    SIM.physicsList = "FTFP_BERT_EMZ"

## FourVector of translation for the Smearing of the Vertex position: x y z t
SIM.vertexOffset = [0.0, 0.0, 0.0, 0.0]
## FourVector of the Sigma for the Smearing of the Vertex position: x y z t
SIM.vertexSigma = [0.0, 0.0, 0.0, 0.0]

## Overwrite Geant4 actions

SIM.action.event = []  # default event action

SIM.action.run = []

SIM.action.stack = []

SIM.action.step = []

SIM.action.track = []

SIM.action.mapActions["DCH_v2"] = "Geant4VoidSensitiveAction"  # drift chamber
SIM.action.calo = "Geant4VoidSensitiveAction"
SIM.action.mapActions["DRcalo"] = "Geant4VoidSensitiveAction"

SIM.action.tracker = (
    "Geant4TrackerWeightedAction",
    {"HitPositionCombination": 2, "CollectSingleDeposits": False},
)

## List of patterns matching sensitive detectors of type Tracker.
SIM.action.trackerSDTypes = ["tracker"]
SIM.filter.tracker = "edep1kev"  #  default for tracking sensitive detectors

##  set a different tracker action for silicon detectors that are thin and require single deposits for accurate delta electron simulation
for region in simulateThinSiliconDetectors:
    SIM.action.mapActions[region] = (
        "Geant4TrackerWeightedAction",
        {"HitPositionCombination": 2, "CollectSingleDeposits": True},
    )
    SIM.filter.mapDetFilter[region] = "edep0"

################################################################################
## Configuration for the output levels of DDG4 components
################################################################################

SIM.output.geometry = 2
SIM.output.inputStage = 3  # input sources
SIM.output.kernel = 3  # Geant4 kernel
SIM.output.part = 3  # ParticleHandler
SIM.output.random = 6  # Eandom Number Generator setup

################################################################################
## Configuration for Output Files.
# (by default, output plugin is determined from output file extension)
################################################################################

## Generator Statuses that are used to mark unstable particles that should decay inside of Geant4.
SIM.physics.alternativeDecayStatuses = set()

################################################################################
## Properties for the random number generator
################################################################################

## If True, calculate random seed for each event basedon eventID and runID
## Allows reproducibility even whenSkippingEvents
SIM.random.enableEventSeed = False
SIM.random.file = None
SIM.random.luxury = 1
SIM.random.replace_gRandom = True
SIM.random.seed = 10
SIM.random.type = None


SIM.ui.commandsConfigure = []
SIM.ui.commandsInitialize = []
SIM.ui.commandsPostRun = []
SIM.ui.commandsPreRun = []
SIM.ui.commandsTerminate = []
