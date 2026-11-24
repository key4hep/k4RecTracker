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
# Test for the TrackMerger algorithm: run the full ILD reconstruction chain, but
# with tracking only, on previously simulated gun muons and verify the merged
# tracks.
#
# Usage: test_track_merger.sh <detector model>
#
# Has to be run from within the ILDConfig/StandardConfig/production directory
# and expects the <detector model>_SIM.edm4hep.root input file to be present
# there.

set -e

if [ -z "$1" ]; then
  echo "Usage: $0 <detector model>" >&2
  exit 1
fi

det_model=$1
script_dir=$(cd "$(dirname "$0")" && pwd)

k4run ILDReconstruction.py \
    --inputFiles=${det_model}_SIM.edm4hep.root \
    --outputFileBase=${det_model} \
    --detectorModel=${det_model} \
    --num-events -1

python3 "${script_dir}/verify_merged_tracks.py" "${det_model}_REC.edm4hep.root"
