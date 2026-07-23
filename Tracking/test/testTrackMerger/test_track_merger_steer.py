"""
Steering file for the TrackMerger test.

Can be run standalone with:

    k4run test_track_merger_steer.py

Input/output file names default to sensible values but are exposed via CLI
so the test driver can override them.
"""

from Configurables import TrackMerger
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser

# ---------------------------------------------------------------------------
# Collection names and matching tolerances.
# Must be kept in sync with the constants in test_track_merger.py.
# ---------------------------------------------------------------------------
INNER_COLL = "InnerTracks"
OUTER_COLL = "OuterTracks"
OUT_COLL = "MyCandidateMergedTracks"

D0_TOL = 0.5
Z0_TOL = 2.5
PHI_TOL = -1.0
OMEGA_TOL = -1.0
TAN_LAMBDA_TOL = -1.0

parser.add_argument("--input", default="input.root", help="Input EDM4hep file")
parser.add_argument("--output", default="output.root", help="Output EDM4hep file")
args = parser.parse_known_args()[0]

iosvc = IOSvc()
iosvc.Input = args.input
iosvc.Output = args.output

merger = TrackMerger(
    "TrackMerger",
    InputInnerTracks=INNER_COLL,
    InputOuterTracks=OUTER_COLL,
    OutTracks=OUT_COLL,
    Greedy=True,
    D0Tolerance=D0_TOL,
    Z0Tolerance=Z0_TOL,
    PhiTolerance=PHI_TOL,
    OmegaTolerance=OMEGA_TOL,
    TanLambdaTolerance=TAN_LAMBDA_TOL,
)

ApplicationMgr(
    TopAlg=[merger],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[iosvc],
)
