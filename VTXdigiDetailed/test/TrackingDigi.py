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
from Gaudi.Configuration import WARNING, DEBUG, INFO
from Configurables import VTXdigitizerDetailed

VTXBarrelDigitizer = VTXdigitizerDetailed(
    "VTXBarrelDigitizer",
    inputSimHits="VertexBarrelCollection",
    outputDigiHits="VXDTrackerHits",
    outputSimDigiAssociation="VXDTrackerHitRelations",
    detectorName="VertexBarrel",
    PixSizePhi=[0.020, 0.020, 0.020, 0.020, 0.020, 0.020],  # in mm
    PixSizeTheta=[0.020, 0.020, 0.020, 0.020, 0.020, 0.020],  # in mm
    tResolution=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    Threshold=100.0,
    ThresholdSmearing=20.0,
    OutputLevel=WARNING,
    DebugHistos=False,
    DebugFileName="Debug_VTXBarrelDigitizer.root",
)

VTXEndcapDigitizer = VTXdigitizerDetailed(
    "VTXEndcapDigitizer",
    inputSimHits="VertexEndcapCollection",
    outputDigiHits="VXDEndcapTrackerHits",
    outputSimDigiAssociation="VXDEndcapTrackerHitRelations",
    detectorName="VertexEndcap",
    PixSizePhi=[0.020, 0.020, 0.020, 0.020, 0.020, 0.020],  # in mm
    PixSizeTheta=[0.020, 0.020, 0.020, 0.020, 0.020, 0.020],  # in mm
    tResolution=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    Threshold=100.0,
    ThresholdSmearing=20.0,
    OutputLevel=WARNING,
    DebugHistos=False,
    DebugFileName="Debug_VTXEndcapDigitizer.root",
)

InnerTrackerBarrelDigitizer = VTXdigitizerDetailed(
    "InnerTrackerBarrelDigitizer",
    inputSimHits="InnerTrackerBarrelCollection",
    outputDigiHits="ITrackerHits",
    outputSimDigiAssociation="InnerTrackerBarrelHitsRelations",
    detectorName="InnerTrackerBarrel",
    PixSizePhi=[0.050, 0.050, 0.050],  # in mm
    PixSizeTheta=[0.300, 0.300, 0.300],  # in mm
    tResolution=[0.0, 0.0, 0.0],
    Threshold=100.0,
    ThresholdSmearing=20.0,
    OutputLevel=WARNING,
    DebugHistos=False,
    DebugFileName="Debug_InnerTrackerBarrelDigitizer.root",
)

InnerTrackerEndcapDigitizer = VTXdigitizerDetailed(
    "InnerTrackerEndcapDigitizer",
    inputSimHits="InnerTrackerEndcapCollection",
    outputDigiHits="ITrackerEndcapHits",
    outputSimDigiAssociation="InnerTrackerEndcapHitsRelations",
    detectorName="InnerTrackerEndcap",
    PixSizePhi=[0.020, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050],  # in mm
    PixSizeTheta=[0.020, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300],  # in mm
    tResolution=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    Threshold=100.0,
    ThresholdSmearing=20.0,
    OutputLevel=WARNING,
    DebugHistos=False,
    DebugFileName="Debug_InnerTrackerEndcapDigitizer.root",
)

OuterTrackerBarrelDigitizer = VTXdigitizerDetailed(
    "OuterTrackerBarrelDigitizer",
    inputSimHits="OuterTrackerBarrelCollection",
    outputDigiHits="OTrackerHits",
    outputSimDigiAssociation="OuterTrackerBarrelHitsRelations",
    detectorName="OuterTrackerBarrel",
    PixSizePhi=[0.050, 0.050, 0.050],  # in mm
    PixSizeTheta=[0.300, 0.300, 0.300],  # in mm
    tResolution=[0.0, 0.0, 0.0],
    Threshold=100.0,
    ThresholdSmearing=20.0,
    OutputLevel=WARNING,
    DebugHistos=False,
    DebugFileName="Debug_OuterTrackerBarrelDigitizer.root",
)

OuterTrackerEndcapDigitizer = VTXdigitizerDetailed(
    "OuterTrackerEndcapDigitizer",
    inputSimHits="OuterTrackerEndcapCollection",
    outputDigiHits="OTrackerEndcapHits",
    outputSimDigiAssociation="OuterTrackerEndcapHitsRelations",
    detectorName="OuterTrackerEndcap",
    PixSizePhi=[0.050, 0.050, 0.050, 0.050],  # in mm
    PixSizeTheta=[0.300, 0.300, 0.300, 0.300],  # in mm
    tResolution=[0.0, 0.0, 0.0, 0.0],
    Threshold=100.0,
    ThresholdSmearing=20.0,
    OutputLevel=WARNING,
    DebugHistos=False,
    DebugFileName="Debug_OuterTrackerEndcapDigitizer.root",
)

TrackingDigiSequence = [
    VTXBarrelDigitizer,
    VTXEndcapDigitizer,
    InnerTrackerBarrelDigitizer,
    InnerTrackerEndcapDigitizer,
    OuterTrackerBarrelDigitizer,
    OuterTrackerEndcapDigitizer,
]
