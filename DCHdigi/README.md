# Drift chamber (DCH) digitizers

## DCHdigi_v01

* Each simulated hit is transformed into a digitized hit. The digitized hit position is the projection of the simulated hit position onto the sense wire (at the center of the cell)
* Smearing of the digitized hit position along the wire and radially is done according to the input parameter values (`zResolution_mm` and `xyResolution_mm`, respectively)
* The digitized hit adds new information: number of clusters and their size, which are derived from precalculated distributions contained in an input file specified by the parameter `fileDataAlg`. The method and distributions corresponds to the option 3 described in F. Cuna et al, arXiv:2105.07064
* It requires that the cellID contain the layer and number of cell within the layer (nphi). It does not matter if the segmentation comes from geometrical segmentation by using twisted tubes and hyperboloids (and the cellID is created out of volume IDs), or the segmentation is virtual DD4hep segmentation
* New digitized hit class is used as an EDM4hep data extension, to be integrated into EDM4hep
* Debug histograms are created if `create_debug_histograms` option is enabled (output file name can be given)
* Stand alone test run simulation of the drift chamber based on twisted tubes, and then apply the digitizer. Dedicated directory with all the files needed is given in `DCHdigi/test/test_DCHdigi/`
* Random number generator uses the seeds calculated on an event basis by the UID service, from the podio header information (run/event number)
* This digitizer is meant to be used with `DriftChamber_o1_v02` and is expected to work for the upcoming `DriftChamber_o1_v03`

## DCHsimpleDigitizer

* Similar to DCHdigi_v01, but it uses the placement matrix of the wires to apply the smearing along/perpendicular to the sense wire
* No cluster counting information is added into the digitized output
* It relies on a dedicated data extension, similar to DCHdigi_v01
