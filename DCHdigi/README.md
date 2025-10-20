# Drift chamber (DCH) digitizers

## DCHdigi_v01

* Each simulated hit is transformed into a digitized hit. The digitized hit position is the projection of the simulated hit position onto the sense wire (at the center of the cell)
* Smearing of the digitized hit position along the wire and radially is done according to the input parameter values (`zResolution_mm` and `xyResolution_mm`, respectively)
* The digitized hit adds dNdx information if flag `calculate_dndx` is enabled (default not). This information consist on number of clusters and their size, which are derived from precalculated distributions contained in an input file specified by the parameter `fileDataAlg`. The method and distributions corresponds to the option 3 described in F. Cuna et al, arXiv:2105.07064
* It requires that the cellID contain the layer and number of cell within the layer (nphi). It does not matter if the segmentation comes from geometrical segmentation by using twisted tubes and hyperboloids (and the cellID is created out of volume IDs), or the segmentation is virtual DD4hep segmentation
* New digitized hit class is used as an EDM4hep data extension, to be integrated into EDM4hep
* Debug histograms are created if `create_debug_histograms` option is enabled (output file name can be given)
* Stand alone test run simulation of the drift chamber based on twisted tubes, and then apply the digitizer. Dedicated directory with all the files needed is given in `DCHdigi/test/test_DCHdigi/`
* Random number generator uses the seeds calculated on an event basis by the UID service, from the podio header information (run/event number)
* This digitizer is meant to be used with `DriftChamber_o1_v02` from k4geo and is expected to work for the upcoming `DriftChamber_o1_v03`

## DCHdigi_v02
* SimHits in one cell are grouped in 'trains' (by time of arrival at the readout), with each train creating one DigiHit (for a single particle traversing one cell, this typically leads to one DigiHit in the cell). Trains are separated by cell dead time (`Deadtime_ns`)
* New digitized hit class is used as an EDM4hep data extension, to be integrated into EDM4hep
* Smearing of the digitized hit position along the wire and radially is done according to the input parameter values (`zResolution_mm` and `xyResolution_mm`, respectively). The digitized hit position is the projection of the simulated hit position onto the sense wire (at the center of the cell)
* Time of DigiHit consists of SimHit creation time +  drift time to wire + signal travelling time along the wire to the readout
* The readout window can be defined via two `Gaudi::Property` members (start time and duration) to filter any DigiHits with associated time outside this window
* dN/dx information is added via Delphes parametrisation. The cluster size is not calculated at the moment, only the total number of clusters for the DigiHit
* It requires that the cellID contain the layer and number of cell within the layer (nphi). It does not matter if the segmentation comes from geometrical segmentation by using twisted tubes and hyperboloids (and the cellID is created out of volume IDs), or the segmentation is virtual DD4hep segmentation
* Random number generator uses the seeds calculated on an event basis by the UID service, from the podio header information (run/event number)
* This digitizer is meant to be used with `DriftChamber_o1_v02` from k4geo and is expected to work for the upcoming `DriftChamber_o1_v03`


## DCHsimpleDigitizerExtendedEdm

* Algorithm for creating digitized drift chamber hits (based on edm4hep::TrackerHit3D) from edm4hep::SimTrackerHit. Resolution along z and xy (distance to the wire) has to be specified. The smearing is applied in the wire reference frame, by means of the placement matrix of the wires
* No cluster counting information is added into the digitized output
* It relies on a dedicated data extension, similar to DCHdigi_v01
* This digitizer is meant to be used with `DriftChamber_o1_v01` from k4geo. Deprecated.
