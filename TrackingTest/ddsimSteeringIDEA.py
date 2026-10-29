# ddsim steering for producing IDEA_o1_v04 simulation input for the tracking
# comparison tests. Derived from the IDEA simulation steering used in the IDEA
# tracking data-production pipelines.
#
# The calorimeters (including the fiber dual-readout calorimeter, which has no
# default sensitive action and would otherwise abort the setup) and the muon
# system are made void to keep the test fast; the drift chamber uses the plain
# tracker action as in the IDEA production steering.
from DDSim.DD4hepSimulation import DD4hepSimulation

SIM = DD4hepSimulation()

# Skip the calorimeters entirely
SIM.action.calo = "Geant4VoidSensitiveAction"
SIM.action.mapActions["DRcalo"] = "Geant4VoidSensitiveAction"

# Set the drift chamber action
SIM.action.mapActions["DCH_v2"] = "Geant4TrackerAction"

# Set the default tracker action
SIM.action.tracker = (
    "Geant4TrackerWeightedAction",
    {"HitPositionCombination": 2, "CollectSingleDeposits": False},
)
SIM.action.trackerSDTypes = ["tracker"]

# Magnetic field stepper configuration used by the IDEA simulation
SIM.field.delta_chord = 0.25
SIM.field.delta_intersection = 0.001
SIM.field.delta_one_step = 0.01
SIM.field.eps_max = 0.001
SIM.field.eps_min = 5e-05
SIM.field.equation = "Mag_UsualEqRhs"
SIM.field.largest_step = 10000.0
SIM.field.min_chord_step = 0.01
SIM.field.stepper = "ClassicalRK4"
