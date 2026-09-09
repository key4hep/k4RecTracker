"""
Steering file for the SimTrackerHitCellMerger test.

Run with:

    k4run test_simtrackerhit_cell_merger_steer.py --input input.root --output output.root

One instance is scheduled per value of the MultipleTrackHandling property, plus one more for the
"Average" RepresentativeKinematics, so that every choice of both properties is exercised in a
single job.

The collection name properties are declared by the k4FWCore KeyValues of the transformer and are
therefore list valued, so they have to be given as lists rather than as bare strings.

Input/output file names default to sensible values but can be overridden via CLI, e.g. by the
CTest setup in CMakeLists.
"""

from Configurables import SimTrackerHitCellMerger
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser

# ---------------------------------------------------------------------------
# Collection names.
# Must be kept in sync with the constants in test_simtrackerhit_cell_merger.py.
# ---------------------------------------------------------------------------
INPUT_COLL = "SimTrackerHits"
OUT_COLL_ALL = "MergedHitsAll"
OUT_COLL_PRIMARY = "MergedHitsPrimaryOnly"
OUT_COLL_PER_TRACK = "MergedHitsPerTrack"
OUT_COLL_SINGLE_TRACK = "MergedHitsSingleTrackCells"
OUT_COLL_AVERAGE = "MergedHitsAllAveraged"

parser.add_argument("--input", default="input.root", help="Input EDM4hep file")
parser.add_argument("--output", default="output.root", help="Output EDM4hep file")
args = parser.parse_known_args()[0]

iosvc = IOSvc()
iosvc.Input = args.input
iosvc.Output = args.output

merger_all = SimTrackerHitCellMerger(
    "SimTrackerHitCellMergerAll",
    InputSimTrackerHits=[INPUT_COLL],
    OutputSimTrackerHits=[OUT_COLL_ALL],
    MultipleTrackHandling="SumAll",
)

merger_primary = SimTrackerHitCellMerger(
    "SimTrackerHitCellMergerPrimaryOnly",
    InputSimTrackerHits=[INPUT_COLL],
    OutputSimTrackerHits=[OUT_COLL_PRIMARY],
    MultipleTrackHandling="PrimaryOnly",
)

merger_per_track = SimTrackerHitCellMerger(
    "SimTrackerHitCellMergerPerTrack",
    InputSimTrackerHits=[INPUT_COLL],
    OutputSimTrackerHits=[OUT_COLL_PER_TRACK],
    MultipleTrackHandling="PerTrack",
)

merger_single_track = SimTrackerHitCellMerger(
    "SimTrackerHitCellMergerSingleTrackCells",
    InputSimTrackerHits=[INPUT_COLL],
    OutputSimTrackerHits=[OUT_COLL_SINGLE_TRACK],
    MultipleTrackHandling="SkipMultiTrackCells",
)

merger_average = SimTrackerHitCellMerger(
    "SimTrackerHitCellMergerAllAveraged",
    InputSimTrackerHits=[INPUT_COLL],
    OutputSimTrackerHits=[OUT_COLL_AVERAGE],
    MultipleTrackHandling="SumAll",
    RepresentativeKinematics="Average",
)

ApplicationMgr(
    TopAlg=[merger_all, merger_primary, merger_per_track, merger_single_track, merger_average],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[iosvc],
)
