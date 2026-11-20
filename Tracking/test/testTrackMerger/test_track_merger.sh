#
# Copyright (c) 2020-2024 Key4hep-Project.
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
# ==============================
#
# Test for the TrackMerger algorithm: generate gun muons, then run the full ILD reconstruction chain, but with tracking only

det_model=ILD_FCCee_v01
det_path=${K4GEO}/FCCee/ILD_FCCee/compact/${det_model}/${det_model}.xml
ddsim --enableGun --gun.distribution "cos(theta)" --gun.energy "5*GeV" --gun.particle "mu-" --gun.multiplicity 10 --numberOfEvents 3 --outputFile ${det_model}_SIM.edm4hep.root --random.enableEventSeed --random.seed 42 --compactFile ${det_path}
k4run ILDReconstruction.py --inputFiles=${det_model}_SIM.edm4hep.root --outputFileBase=${det_model} --detectorModel=${det_model} --trackingOnly
python3 verify_merged_tracks.py "${det_model}_REC.edm4hep.root"

det_model=ILD_FCCee_v02
det_path=${K4GEO}/FCCee/ILD_FCCee/compact/${det_model}/${det_model}.xml
ddsim --enableGun --gun.distribution "cos(theta)" --gun.energy "5*GeV" --gun.particle "mu-" --gun.multiplicity 10 --numberOfEvents 3 --outputFile ${det_model}_SIM.edm4hep.root --random.enableEventSeed --random.seed 42 --compactFile ${det_path}
k4run ILDReconstruction.py --inputFiles=${det_model}_SIM.edm4hep.root --outputFileBase=${det_model} --detectorModel=${det_model} --trackingOnly
python3 verify_merged_tracks.py "${det_model}_REC.edm4hep.root"
