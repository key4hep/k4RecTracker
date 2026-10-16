"""
Steering file for the TrackMerger test.

Run with:

    k4run test_track_merger_steer.py --mode absolute --input input.root --output output.root

--mode selects between absolute differences (default) and significances. The
significance mode loads field_only.xml through the GeoSvc to provide the
magnetic field needed for the pT significance.

Input/output file names default to sensible values but can be overridden
via CLI, e.g. by the CTest setup in CMakeLists.
"""

import os

from Configurables import GeoSvc, TrackMerger
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser

# ---------------------------------------------------------------------------
# Collection names and matching tolerances.
# Must be kept in sync with the constants in test_track_merger.py.
# ---------------------------------------------------------------------------
INNER_COLL = "InnerTracks"
OUTER_COLL = "OuterTracks"
OUT_COLL = "MyCandidateMergedTracks"

# Absolute mode
D0_TOL = 0.5
Z0_TOL = 2.5
PHI_TOL = -1.0
OMEGA_TOL = -1.0
TAN_LAMBDA_TOL = -1.0

# Significance mode
D0_SIG = 3.0
Z0_SIG = 3.0
PHI_SIG = 3.0
PT_SIG = 2.0
TAN_LAMBDA_SIG = 3.0

parser.add_argument(
    "--mode", choices=["absolute", "significance"], default="absolute", help="Matching mode"
)
parser.add_argument("--input", default="input.root", help="Input EDM4hep file")
parser.add_argument("--output", default="output.root", help="Output EDM4hep file")
args = parser.parse_known_args()[0]

iosvc = IOSvc()
iosvc.Input = args.input
iosvc.Output = args.output

ext_svc = [iosvc]

if args.mode == "absolute":
    cuts = dict(
        UseSignificance=False,
        D0Tolerance=D0_TOL,
        Z0Tolerance=Z0_TOL,
        PhiTolerance=PHI_TOL,
        OmegaTolerance=OMEGA_TOL,
        TanLambdaTolerance=TAN_LAMBDA_TOL,
    )
else:
    geoservice = GeoSvc("GeoSvc")
    test_dir = os.path.dirname(os.path.abspath(__file__))
    geoservice.detectors = [os.path.join(test_dir, "field_only.xml")]
    geoservice.EnableGeant4Geo = False
    ext_svc.insert(0, geoservice)
    cuts = dict(
        UseSignificance=True,
        D0Significance=D0_SIG,
        Z0Significance=Z0_SIG,
        PhiSignificance=PHI_SIG,
        PtSignificance=PT_SIG,
        TanLambdaSignificance=TAN_LAMBDA_SIG,
    )

merger = TrackMerger(
    "TrackMerger",
    InputInnerTracks=INNER_COLL,
    InputOuterTracks=OUTER_COLL,
    OutTracks=OUT_COLL,
    Greedy=True,
    **cuts,
)

ApplicationMgr(
    TopAlg=[merger],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=ext_svc,
)
