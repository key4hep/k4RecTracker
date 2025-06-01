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