#
# Copyright (c) 2014-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from Gaudi.Configuration import WARNING
from Configurables import MarlinProcessorWrapper


OverlayParameters = {
    "MCParticleCollectionName": ["MCParticle"],
    "MCPhysicsParticleCollectionName": ["MCPhysicsParticles"],
    "Delta_t": ["20"],
    "NBunchtrain": ["20"],
    "Collection_IntegrationTimes": [
        "VertexBarrelCollection", "380",
        "VertexEndcapCollection", "380",
        "InnerTrackerBarrelCollection", "380",
        "InnerTrackerEndcapCollection", "380",
        "OuterTrackerBarrelCollection", "380",
        "OuterTrackerEndcapCollection", "380",
        "ECalBarrelCollection", "380",
        "ECalEndcapCollection", "380",
        "HCalBarrelCollection", "380",
        "HCalEndcapCollection", "380",
        "HCalRingCollection", "380",
        "YokeBarrelCollection", "380",
        "YokeEndcapCollection", "380",
        "LumiCalCollection", "380",
        "ECalBarrelEta", "380",  # for CLD_o4_v05 with LAr
        "ArcCollection", "380",  # for CLD_o3_v01 with ARC
     ],
    "PhysicsBX": ["1"],
    "Poisson_random_NOverlay": ["false"],
    "RandomBx": ["false"],
    "TPCDriftvelocity": ["0.05"],
    # TODO: add argument, also add argument to turn overlay on
    "BackgroundFileNames": ["pairs_Z_sim.slcio"],
}
Overlay = {}

Overlay["False"] = MarlinProcessorWrapper("OverlayFalse")
Overlay["False"].OutputLevel = WARNING
Overlay["False"].ProcessorType = "OverlayTimingGeneric"
Overlay["False"].Parameters = OverlayParameters.copy()
Overlay["False"].Parameters |= {
                           "BackgroundFileNames": [],
                           "NBunchtrain": ["0"],
                           "NumberBackground": ["0."],
                           }

# XXX: Caution, this probably needs an update
Overlay["91GeV"] = MarlinProcessorWrapper("Overlay91GeV")
Overlay["91GeV"].OutputLevel = WARNING
Overlay["91GeV"].ProcessorType = "OverlayTimingGeneric"
Overlay["91GeV"].Parameters = OverlayParameters.copy()
Overlay["91GeV"].Parameters |= {
                           "NumberBackground": ["1."],
                           }

# XXX: Caution, this probably needs an update
Overlay["365GeV"] = MarlinProcessorWrapper("Overlay365GeV")
Overlay["365GeV"].OutputLevel = WARNING
Overlay["365GeV"].ProcessorType = "OverlayTimingGeneric"
Overlay["365GeV"].Parameters = OverlayParameters.copy()
Overlay["365GeV"].Parameters |= {
                            "Delta_t": ["3396"],
                            "NBunchtrain": ["3"],
                            "NumberBackground": ["1."],
                            }

OverlaySequence = [
    Overlay[CONFIG["Overlay"]],
]
