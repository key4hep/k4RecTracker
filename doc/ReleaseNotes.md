# v00.08.00

* 2026-10-16 Victor Schwan ([PR#113](https://github.com/Key4hep/k4RecTracker/pull/113))
  - Introduce significance as matching criterion for the `TrackMerger`
  - Update matching defaults to match defaults from `ClonesAndSplitTracksFinder` in `CLDConfig`
  - Change default behavior to use the significance as matching criterion

* 2026-10-15 Giovanni Marchiori ([PR#117](https://github.com/Key4hep/k4RecTracker/pull/117))
  - Retrieve sensor thickness from volume rather than from ZPlanarData and ZPetalData extensions.

* 2026-10-15 andread3vita ([PR#116](https://github.com/Key4hep/k4RecTracker/pull/116))
  - Add EmbeddingCoordinates to the Gaudi properties of GGTFTrackFinder to support models trained with different embedding space dimensions.
  - Support extra features in the output of the model. The extra features are dropped and not used by the clustering, but the inference does not crash.

* 2026-10-14 andread3vita ([PR#112](https://github.com/Key4hep/k4RecTracker/pull/112))
  - Reset the global TGeoManager navigator to a canonical state before each fit. The navigator state is process-global and is modified by material lookups during track extrapolation. Without this reset, the material assigned to points located on geometry boundaries may depend on the previously fitted track, potentially from an earlier event, making the fit results dependent on processing order.
  - Remove isPositiveSemiDefinite from utils.

* 2026-10-07 Juan Miguel Carceller ([PR#109](https://github.com/Key4hep/k4RecTracker/pull/109))
  - TracksFromGenParticles: fix a bug where, with `ExtrapolateToECal=False`, the ECAL geometry members were left uninitialized but still read for every track, so tracks could randomly get calorimeter track states extrapolated with garbage geometry. All the members are now default-initialized and the extrapolation is only done when `ExtrapolateToECal` is set.
  - TracksFromGenParticles: return `StatusCode::FAILURE` with an error message instead of throwing when the calorimeter geometry lookup fails, remove the redundant `retrieved` flag and use separate pointers for the DRC, barrel and endcap extensions.
  - TracksFromGenParticles: store the system encoder as a member instead of a raw pointer, use a 64-bit type for the cellID (it was truncated to `uint`) and add missing includes.

* 2026-10-06 Juan Miguel Carceller ([PR#110](https://github.com/Key4hep/k4RecTracker/pull/110))
  - Fix uninitialized variable warnings in VTXdigitizerDetailed

* 2026-10-06 andread3vita ([PR#102](https://github.com/Key4hep/k4RecTracker/pull/102))
  - Refactored GGTFTrackFinder by splitting the reconstruction workflow into dedicated functions for input preparation, ONNX inference, and track building, improving readability and maintainability.
  - Introduced a HitBatch structure to consistently store flattened model inputs together with the corresponding hit type, collection index, and hit index.
  - Fixed the bookkeeping of wire-hit collection indices, which previously caused hits from multiple wire input collections to be associated with the wrong collection.
  - Added filtering of noise/unclustered hits: cluster IDs less than or equal to zero no longer produce output tracks, and only successfully reconstructed tracks are written to the output collection.
  - Removed the use of track type 1 when saving reconstructed tracks.
  - The number of hits is now checked before reserving and filling the model input buffers, avoiding unnecessary memory allocation for events that will be rejected because numHits > kMaxHits.
  - Improved memory management of ONNX input and output names by storing owned copies instead of releasing the memory returned by ONNX Runtime, fixing the previous memory leak.
  - Reworked wire-hit geometry construction and input feature preparation with explicit bookkeeping of the original source collections and hit indices.
  - Made the processed-event counter thread-safe by replacing the mutable integer counter with an atomic counter.
  - Removed unused standard-library, Gaudi, k4Interface, EDM4HEP, ATen, and ONNX Runtime includes.
  - Simplified ONNX Runtime resource ownership using RAII-managed environment, session, memory information, and owned input/output names.
  - Added ONNX model and inference validation, including checks for a valid model path, exactly one input and one output, consistent input buffer and tensor dimensions, and valid output tensor size, with improved error handling during ONNX Runtime initialization and inference.

* 2026-10-05 Jona Dilg ([PR#111](https://github.com/Key4hep/k4RecTracker/pull/111))
  - Fix typo in CMakeLists of ARCdigi and VTXdigi that would install header files into `include/@{CMAKE_PROJECT_NAME}` instead of expanding `${CMAKE_PROJECT_NAME}`

* 2026-10-01 Juan Miguel Carceller ([PR#106](https://github.com/Key4hep/k4RecTracker/pull/106))
  - Fix deprecated getTrackState usage in TrackMerger component by using `getTrackState` where it is needed.

* 2026-09-30 Juan Miguel Carceller ([PR#105](https://github.com/Key4hep/k4RecTracker/pull/105))
  - Remove the DD4hep rpath and compiler flags setting from CMakeLists.txt. These are set by Key4hepConfig.cmake; the DD4hep setting is older than the inclusion of Key4hepConfig.cmake.
   - Remove a duplicated `CMAKE_MODULE_PATH` setting.
   - Require onnxruntime 1.16.0, the first version that installs a CMake config with the `onnxruntime::onnxruntime` target.

* 2026-09-30 andread3vita ([PR#101](https://github.com/Key4hep/k4RecTracker/pull/101))
  - Fix the helix-to-Cartesian conversion and track initialization:
    - Correct the units used in the Jacobian to ensure consistency between GENFIT and EDM4hep.
    - Correct the unit conversion of z0.
    - Prevent a singular initial covariance in z0 when the track seed is close to the reference point.
  - Improve the construction of EDM4hep track states:
    - Use GENFIT native extrapolation methods to propagate fitted states to the point of closest approach to the requested reference point.
    - Use the first and last digitized hits as reference points when constructing the AtFirstHit and AtLastHit track states.
    - Make the AtIP extrapolation consistent for prompt tracks.
    - Document that the AtIP track state is not meaningful for displaced tracks when the displaced-track configuration uses a non-IP reference point.
  - Make SetVPPosition part of the public GenfitTrack API, allowing the vertex/reference position used for track-state construction to be configured.
  - Fix hit-to-track associations affected by hit filtering and change the handling of FittedHits:
    - Modify the provided FittedHits collection instead of creating it inside ProcessTrack.
    - Store only successfully fitted/accepted hits; rejected-hit placeholders are no longer inserted, so FittedHits is no longer index-aligned with the input hits.
    - Fill fitted hits only when the GENFIT fit and the required track-state propagation complete successfully.
  - Improve failure handling after track fitting: failures while retrieving or propagating a fitted state are caught per track and no longer terminate processing of the full event.
  - Manage m_genfitTrack with std::unique_ptr, ensuring automatic cleanup and removing manual lifetime management.

* 2026-09-28 Federico Meloni ([PR#97](https://github.com/Key4hep/k4RecTracker/pull/97))
  - Add K4RECTRACKER_TRACKING CMake flag (default=ON) gating the compilation of the FCCee tracking

* 2026-09-23 Andreas Loeschcke Centeno ([PR#103](https://github.com/Key4hep/k4RecTracker/pull/103))
  - Fixed the `DCHdigi_v02` alias for the `WireTrackerDigi_v01`

* 2026-09-15 Andreas Loeschcke Centeno ([PR#92](https://github.com/Key4hep/k4RecTracker/pull/92))
  - Adapted and changed logic of accepted $\beta\gamma$ range in `TrackdNdxDelphesBased`

* 2026-09-06 Juan Miguel Carceller ([PR#99](https://github.com/Key4hep/k4RecTracker/pull/99))
  - Remove goto statements from TrackdNdxDelphesBased.

* 2026-09-03 Juan Miguel Carceller ([PR#98](https://github.com/Key4hep/k4RecTracker/pull/98))
  - Adapt TrackdNdxDelphesBased to the new EDM4hep TrackState lookup API.

* 2026-08-26 Juan Miguel Carceller ([PR#96](https://github.com/Key4hep/k4RecTracker/pull/96))
  - Link Tracking against the imported Pandora SDK target when available for compatibility with Pandora v5

* 2026-08-21 SanghyunKo ([PR#95](https://github.com/Key4hep/k4RecTracker/pull/95))
  - Add a printout message when TracksFromGenParticles failed to find the DR calo geometry

* 2026-08-07 David ([PR#93](https://github.com/Key4hep/k4RecTracker/pull/93))
  - Change geometry returned by geosvc to be used as const in GenfitTrackFitter.cpp

* 2026-08-03 Victor Schwan ([PR#87](https://github.com/Key4hep/k4RecTracker/pull/87))
  - Adds a new `TrackMerger` Gaudi transformer that merges tracks from two input collections based on compatible track-state parameters at their adjoining hit
  - While generic and applicable to merging any two track collections, it was designed to combine silicon-tracker tracks with TPC tracks reconstructed by `Clupatra` for ILD
  - The exact matching criterion in track-parameter space is not yet optimized and may need tuning
  - Includes self-contained CTest `test_TrackMerger*` that generates synthetic matching/non-matching track pairs and verifies correct merging behavior

* 2026-07-30 Juan Miguel Carceller ([PR#89](https://github.com/Key4hep/k4RecTracker/pull/89))
  - Add all the plugin build directories to `LD_LIBRARY_PATH` and `GAUDI_PLUGIN_PATH` when setting up the test environment, so that plugins are found in a clean build tree to avoid issues with the LCG stacks (where there is always a `k4RecTracker` from the stack.

* 2026-07-27 SanghyunKo ([PR#86](https://github.com/Key4hep/k4RecTracker/pull/86))
  - Fix TracksFromGenParticles - build tracks only for particles that reach the calorimeter

* 2026-06-25 Stefano Franchellucci ([PR#83](https://github.com/Key4hep/k4RecTracker/pull/83))
  - Update code to use WireTracker info instead of DCH, allowing for easier compatibility with other trackers

* 2026-05-09 Juan Miguel Carceller ([PR#82](https://github.com/Key4hep/k4RecTracker/pull/82))
  - Remove the deprecated k4DataSvc and PodioOutput, following https://github.com/key4hep/k4FWCore/pull/392

* 2026-05-05 andread3vita ([PR#79](https://github.com/Key4hep/k4RecTracker/pull/79))
  - Correct the data type of phi0 when saving it in PCAInfo
  - Correct the sign of omega by properly determining the particle charge during initialization
  - Add covariance matrix to fitted trackStates
  - Fix how filtered tracks are saved in the output file: now for each fitted track there is the corresponding filtered track
  - Expose EDM4Hep units to the user

* 2026-04-24 andread3vita ([PR#77](https://github.com/Key4hep/k4RecTracker/pull/77))
  - Rename the track finder functional to align with the current coding style
  - Add GenFit interfaces to enable running the track fitter
  - Implement the GenFit-based track fitter
  - Implement PerfectTrackFinder
  - Add test for the track fitter

* 2026-02-24 Mateusz Jakub Fila ([PR#76](https://github.com/Key4hep/k4RecTracker/pull/76))
  - Remove histogram interface compatibility with Gaudi versions older that v39

# v00.07.00

* 2025-11-11 Juan Miguel Carceller ([PR#73](https://github.com/key4hep/k4RecTracker/pull/73))
  - Avoid warnings about possible dangling references when storing the transformation matrix, similarly to how it was done in https://github.com/HEP-FCC/k4RecCalorimeter/pull/184#issue-3517931900

* 2025-11-05 Thomas Madlener ([PR#71](https://github.com/key4hep/k4RecTracker/pull/71))
  - Add pre-commit config and enforce formatting rules via CI

* 2025-11-05 Andreas Loeschcke Centeno ([PR#67](https://github.com/key4hep/k4RecTracker/pull/67))
  - `DCHdigi_v02`
     - combines SimHits in a cell to one DigiHit (potentially multiple DigiHits, if the SimHits are significantly separated by time)
     - Based on exchange with experts, time and position measurements of DigiHit are governed by first hit arriving at the readout 
     - dN/dx calculation based on delphes parametrisation

* 2025-11-04 Sanghyun Ko ([PR#62](https://github.com/key4hep/k4RecTracker/pull/62))
  - Case handling to enable the track (from GenPtc) extrapolation for IDEA

# v00.06.00

* 2025-10-19 Andrea De Vita ([PR#65](https://github.com/key4hep/k4RecTracker/pull/65))
  - A new functional is available in `Tracking/components/GGTF_tracking.cpp`. It supports track finding using CLD, IDEA, and ALLEGRO detectors. The function returns `extension::Tracks` to ensure compatibility with IDEA and ALLEGRO, whose wire hits are not included in `edm4hep`. To run the track finder, you must specify the path to the trained model via the `SteeringFile`. A simple pre-trained model using IDEA v3 o1 is provided in the `test` folder as an example.
  - The utility functions previously used exclusively in Tracking/components/TracksFromGenParticles.cpp have been moved to Tracking/include/utils.hpp, making them accessible to other functional components as well.
  - `test_TracksFromGenParticles` now includes both the data production and the functional test. These two steps have been integrated into a single test, simplifying the workflow and reducing redundancy.

* 2025-10-15 sss ([PR#70](https://github.com/key4hep/k4RecTracker/pull/70))
  - Fix cmake developer warnings.

* 2025-10-15 sss ([PR#69](https://github.com/key4hep/k4RecTracker/pull/69))
  - Fix an unused variable warning seen with clang.

* 2025-07-28 jmcarcell ([PR#66](https://github.com/key4hep/k4RecTracker/pull/66))
  - Use `REQUIRED` with `find_package` to fail immediately, also avoid duplicated `find_package` calls
  - Link to `k4Interface` and `DDRec`, which were used but not linked to. It was working fine when using the stack but locally I was getting errors

* 2025-07-07 Andreas Loeschcke Centeno ([PR#61](https://github.com/key4hep/k4RecTracker/pull/61))
  - New funcitonal _TrackdNdxDelphesBased_ to attach the dN/dx information to tracks, using the delphes approach:
    - Uses track length and particle βγ to calculate number of clusters for track (calling delphes functions directly)

* 2025-06-20 Andrea De Vita ([PR#63](https://github.com/key4hep/k4RecTracker/pull/63))
  - [DCHdigi_v01] PositionAlongWireError and DistanceToWireError are now populated with the corresponding resolution values, rather than the shift applied to the simulated hit.

* 2025-05-30 jmcarcell ([PR#60](https://github.com/key4hep/k4RecTracker/pull/60))
  - Move DataHandle and MetaDataHandle to the k4FWCore namespace, done in k4FWCore in https://github.com/key4hep/k4FWCore/pull/317
  - Bump the required version of k4FWCore

# v00.05.00

* 2025-05-27 jmcarcell ([PR#52](https://github.com/key4hep/k4RecTracker/pull/52))
  - Limit scope of random distributions to the call operator to prevent race conditions

* 2025-04-15 Andrea De Vita ([PR#58](https://github.com/key4hep/k4RecTracker/pull/58))
  - Corrected the unit in the comment describing spatial resolution (`VTXdigi/include/VTXdigitizer.h`): previously stated as µm, but the correct unit is mm.

* 2025-04-01 Giovanni Marchiori ([PR#57](https://github.com/key4hep/k4RecTracker/pull/57))
  - Print out full particle details only in verbose mode, add succinct summary in debug mode

* 2025-03-24 Giovanni Marchiori ([PR#56](https://github.com/key4hep/k4RecTracker/pull/56))
  - add number of hits per tracking sub detector to tracks from gen-level particles
  - ignore hits from secondaries, that screwed up the calculation of number of hits, of track state at last hit and extrapolation to calorimeter
  - use last hit information for extrapolation to calorimeter
  - remove duplicate algorithms

* 2025-03-24 jmcarcell ([PR#54](https://github.com/key4hep/k4RecTracker/pull/54))
  - Fix failing tests in CI by adding a test that creates an input file

* 2025-03-17 jmcarcell ([PR#53](https://github.com/key4hep/k4RecTracker/pull/53))
  - Use the properties `Input` and `Output` with `IOSvc` instead of the deprecated `input` and `output`

* 2025-02-28 jmcarcell ([PR#51](https://github.com/key4hep/k4RecTracker/pull/51))
  - Set test dependencies for tracking tests that need a simulated file

* 2025-02-28 Alvaro Tolosa Delgado ([PR#50](https://github.com/key4hep/k4RecTracker/pull/50))
  - DCH digi output now contains correct values for DistanceToWireError and PositionAlongWireError. Before it was zero. Fixes https://github.com/key4hep/k4RecTracker/issues/48

* 2025-02-28 Mateusz Jakub Fila ([PR#49](https://github.com/key4hep/k4RecTracker/pull/49))
  - Use `getUniqueID` directly with `EventHeaderCollection` to avoid truncating casts
  - Add seeding of TRandom3 rng in DCHdigi_v01. Make TF1 draw random numbers from instad of global rng.

* 2025-02-28 Giovanni Marchiori ([PR#45](https://github.com/key4hep/k4RecTracker/pull/45))
  - add retrieval of B field and extrapolation to calorimeter to functional that produces tracks from gen particles

* 2025-02-27 jmcarcell ([PR#47](https://github.com/key4hep/k4RecTracker/pull/47))
  - Add a requirement on the version of podio needed to make building on top of
    the latest release fail faster. https://github.com/AIDASoft/podio/pull/714 is
    needed and is only avaliable in 1.2 onwards
  - Add LANGUAGES CXX to CMakeLists.txt to disable checks for the C compiler

# v00.04.00

* 2025-01-31 Giovanni Marchiori ([PR#43](https://github.com/key4hep/k4RecTracker/pull/43))
  - Add a Gaudi::Algorithm to create tracks from generator-level particles with charge!=0 and number of simtrackerhits>0 (temporary work around needed because one can not use the edm4hep <--> lcio converters with IOSvc and functionals)
  - Track states at IP, first tracker hit, last tracker hit and extrapolation to calorimeter are calculated and saved by the algorithm

* 2024-12-29 jmcarcell ([PR#40](https://github.com/key4hep/k4RecTracker/pull/40))
  - Remove the check for TrackerHit3D from edm4hep

* 2024-12-26 jmcarcell ([PR#29](https://github.com/key4hep/k4RecTracker/pull/29))
  - Add REQUIRED for MarlinUtil in CMakeLists.txt

* 2024-12-19 Alvaro Tolosa Delgado ([PR#42](https://github.com/key4hep/k4RecTracker/pull/42))
  - DCHdigi_v01 now uses evolved data extension `SenseWireHit`

* 2024-12-18 BrieucF ([PR#41](https://github.com/key4hep/k4RecTracker/pull/41))
  - Add SenseWireHit data type in extension to validate the content of https://github.com/key4hep/EDM4hep/pull/385

* 2024-11-27 Alvaro Tolosa Delgado ([PR#27](https://github.com/key4hep/k4RecTracker/pull/27))
  - New digitizer for DCH v2 (based on twisted tubes). It calculates cluster information and smears the position along and perpendicular the wire separately. New custom EDM4hep data extension was added to include these data.

* 2024-11-25 Archil Durglishvili ([PR#39](https://github.com/key4hep/k4RecTracker/pull/39))
  - Only a charged genParticle is considered to form a track in TracksFromGenParticles algorithm

* 2024-11-08 jmcarcell ([PR#37](https://github.com/key4hep/k4RecTracker/pull/37))
  - Fix service retrieval after deprecations in Gaudi v39r1, see https://gitlab.cern.ch/gaudi/Gaudi/-/merge_requests/1637

# v00-03-00

* 2024-09-28 jmcarcell ([PR#35](https://github.com/key4hep/k4RecTracker/pull/35))
  - Fix a few warnings after the change Association -> Link

* 2024-09-27 tmadlener ([PR#34](https://github.com/key4hep/k4RecTracker/pull/34))
  - Make sure that histograms keep working / compiling with Gaudi v39

* 2024-09-24 armin-ilg ([PR#33](https://github.com/key4hep/k4RecTracker/pull/33))
  - Update to IDEA_o1_v03 vertex
  - Specify spatial and time resolutions per layer

* 2024-09-10 jmcarcell ([PR#32](https://github.com/key4hep/k4RecTracker/pull/32))
  - Use the Key4hepConfig flag to set the standard, compiler flags and rpath magic.

* 2024-09-07 jmcarcell ([PR#30](https://github.com/key4hep/k4RecTracker/pull/30))
  - Remove the old key4hep nightlies workflow

* 2024-08-08 BrieucF ([PR#28](https://github.com/key4hep/k4RecTracker/pull/28))
  - Add a link between SIM and Digitized VTX hits
  - Change detector model to IDEA_o1_v03 for the VTX tests

* 2024-08-08 BrieucF ([PR#21](https://github.com/key4hep/k4RecTracker/pull/21))
  - Add a transformer to create Tracks out of generated particles
  - Add a consumer that plots the distance between simTrackerHits and a Track

* 2024-07-27 jmcarcell ([PR#26](https://github.com/key4hep/k4RecTracker/pull/26))
  - Add missing k4FWCore::k4Interface when linking, it's needed (although it works fine in the stack) because some files from k4Interface are being included
  - Fix warning about an unused variable

* 2024-07-27 jmcarcell ([PR#25](https://github.com/key4hep/k4RecTracker/pull/25))
  - Change headers and add EventContext in algorithms not to use GaudiAlg
  - Replace `GaudiTool` with `AlgTool`

* 2024-07-04 tmadlener ([PR#23](https://github.com/key4hep/k4RecTracker/pull/23))
  - Remove setting a dummy dE/dx value, since the edm4hep::Track will no longer have it after [EDM4hep#311](https://github.com/key4hep/EDM4hep/pull/311)

* 2024-03-19 jmcarcell ([PR#20](https://github.com/key4hep/k4RecTracker/pull/20))
  - Add missing GaudiAlgLib, which seems to be giving some problems now in the nightlies

* 2024-02-26 BrieucF ([PR#19](https://github.com/key4hep/k4RecTracker/pull/19))
  - Drift chamber digitized hit are now in global coordinate by default
  - Added an edm4hep association between drift chamber digi and simTrackerHit 
  - Add a debug mode in the drift chamber digitizer which produces additional distributions and a digi collection in local coordinates

* 2024-02-23 jmcarcell ([PR#18](https://github.com/key4hep/k4RecTracker/pull/18))
  - Change clone paths in the readme with the new location of this repository

* 2024-02-23 jmcarcell ([PR#17](https://github.com/key4hep/k4RecTracker/pull/17))
  - Use TrackerHit3D when avaliable, needed after https://github.com/key4hep/EDM4hep/pull/252

# v00-02-00

* 2024-02-12 tmadlener ([PR#16](https://github.com/key4hep/k4RecTracker/pull/16))
  - Add a schema_version to the data model extension to keep it working after [podio#556](https://github.com/AIDASoft/podio/pull/556)

# v00-01

* This file is also automatically populated by the tagging script