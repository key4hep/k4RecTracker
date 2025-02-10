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