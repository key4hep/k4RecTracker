# TrackingTest

Comparison of two track-reconstruction pipelines for IDEA_o1_v04 on the same
simulated samples:

- **Pipeline A: GGTF + Genfit.** `GGTFTrackFinder`
  (`/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/Tracking/components/GGTFTrackFinder.cpp`)
  with the ONNX model
  `/afs/cern.ch/work/a/adevita/public/trainingFolder/GGTF_trackFinder_training/conversionModel/smallModelForTest.onnx`,
  followed by `GenfitTrackFitter`
  (`/afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/Tracking/components/GenfitTrackFitter.cpp`).
- **Pipeline B: ConformalTracking + Refit.** `ConformalTracking` in k4Reco using
  only the vertex and silicon wrapper digitized hits, followed by `RefitFinal`
  (`/afs/cern.ch/work/a/adevita/public/workDir/k4Reco/k4Reco/Tracking/src/RefitFinal.cpp`).

Both pipelines read the same digitized file, so the input smearing is identical.

Two samples are used:

- **muon gun**: 10 events, 3 mu- @ 10 GeV from a uniform particle gun
  (`run_all.sh`).
- **Z -> uds at 91 GeV**: Pythia generation + IDEA ddsim with crossing angle
  and beam-spot smearing, following the production pipeline
  `runSequence.sh` (`run_zuds.sh`, currently 1 event).

## Files

| File | Purpose |
|---|---|
| `ddsimSteeringIDEA.py` | ddsim steering for the muon-gun sim (void calorimeters/DRcalo, DCH tracker action, field config) |
| `runIDEA_digi.py` | sim -> digitized hits (`VTXBDigis`, `VTXDDigis`, `SiWrBDigis`, `SiWrDDigis`, `DCH_DigiCollection`) |
| `runIDEA_GGTF_Genfit.py` | digitized -> GGTF track finding + Genfit refit (pipeline A) |
| `runIDEA_ConformalTracking_Refit.py` | digitized -> conformal tracking (VTX+SiWr only) + RefitFinal (pipeline B) |
| `compare_tracking.py` | compares two track outputs against generator truth (`--truth gun` or `--truth stable`) |
| `plot_tracking.py` | residual / hit-content / efficiency plots |
| `run_all.sh` | muon-gun chain (sim -> digi -> A, B -> comparison) |
| `Zcard.cmd`, `pythia.py`, `SteeringFile_IDEA_o1_v04.py` | Z -> uds production-style generation and simulation steering (copied from the production pipeline) |
| `run_zuds.sh` | Z -> uds chain (pythia -> ddsim -> digi -> A, B -> comparison) |
| `setup_env.sh` | runtime environment: key4hep stack + local k4RecTracker and k4Reco builds |

## Running

A key4hep environment must be set up first (the scripts assume `K4GEO` is
defined), then:

```bash
cd /afs/cern.ch/work/a/adevita/public/workDir/k4RecTracker/TrackingTest
bash run_all.sh     # muon gun comparison, ends with the comparison printout
bash run_zuds.sh    # Z -> uds comparison (1 event)
```

## Notes

- The GGTF finder consumes all four silicon hit collections plus the drift
  chamber wire hits (as given in the reference configuration). The conformal
  tracking side uses only vertex + silicon wrapper hits; the hit-content
  difference is reported in the comparison output.
- `compare_tracking.py` converts the reconstructed `omega` into pT assuming
  the nominal IDEA 2 T solenoid field (validated by the agreement with the MC
  pT for matched tracks).
- Matching: track state at the IP within |dphi| < 0.02 rad and
  |d(tanLambda)| < 0.02 of a generator particle (greedy best match). Results
  are stable against widening the window to 100 mrad (see sweep in the
  muon-gun discussion).
