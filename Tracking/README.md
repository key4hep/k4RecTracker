# Key4Hep Tracking Algorithms for FCC-ee detectors

This subfolder contains the implementation of several Tracking tools for FCC-ee detectors using the key4hep framework:

* GGTFTrackFinder
* PlotTrackHitResiduals
* TrackdNdxDelphesBased
* TracksFromGenParticles

## Geometric Graph Track Finding

The **Geometric Graph Track Finding (GGTF)** method is an end-to-end, detector-agnostic approach to tracking pattern recognition. It provides a generalized geometric strategy for track finding that:

1. accommodates multiple sub-detectors with heterogeneous input geometries and tracking technologies,
2. does not require detailed knowledge of the detector geometry or material composition, and
3. avoids reliance on analytical parametrizations of particle trajectories.

In our end-to-end pipeline, hits from all tracking components are directly processed to produce a set of reconstructed tracks. A key innovation is the use of a *geometric algebra representation* of the data, which enables the integration of diverse geometric types. This is combined with a graph neural network, **GATr**, designed to exploit detector symmetries through equivariance.

### Technical Implementation

We implemented a **Geometric Graph Track Finder (GGTF)** in `k4RecTracker/Tracking/components/GGTFTrackFinder.cpp`.
Its workflow can be summarized as follows:

1. **Input Extraction**
   - From the hit collections, a 7-dimensional tensor is built.
   - The components are:
     - **[0–2]**: 3D position of the silicon detector hits (e.g., vertex and wrapper).
     - **[3]**: hit type (`0` = silicon detector hit, `1` = drift chamber hit).
     - **[4–6]**: vector components pointing from the left to the right positions along the circles that identify the drift chamber hits (set to `0` for silicon hits).

2. **Machine Learning Step**
   - The 7-dimensional inputs are mapped into a collection of 4-dimensional points in an embedding space.
   - Each 4D point consists of **3 geometric coordinates** and **1 charge-like component**.
   - Intuitively, this charge can be seen as a potential that **attracts hits of the same cluster** and **repels unrelated ones**.
   - This step is implemented with an [`Ort::Session`](https://onnx.ai/) initialized with a `.onnx` model.

3. **Clustering Step**
   - The 4D embedded points are clustered into sub-collections.
   - Each cluster corresponds **one-to-one** to a reconstructed track.

4. **Track Creation**
   - The identified clusters are converted into tracks.
   - The final results are stored in the collection `OutputTracksGGTF`.

### Dependencies

* ROOT
* PODIO
* Gaudi
* k4FWCore
* ONNX

### Installation

This algorithm can be installed by compiling the whole k4RecTracker project.

### Execution

Run the tagger by including the **GGTF_tracking** algorithm in a steering file like [runTestTrackFinder.py](test/testTrackFinder/runTestTrackFinder.py) and run it like this:

```bash
k4run test/testTrackFinder/runTestTrackFinder.py --modelPath <MODEL_PATH> --tbeta <TBETA> --td <TD>
```

This will return your edm4hep output file with the added `OutputTracksGGTF` collection.

### Retraining a model

When the project is compiled, the latest trained model for **IDEA_o1_v2** is automatically downloaded.
This model has been trained on **Z → qq** events at 91 GeV without background.

It is recommended to re-train the model if:

* you plan to apply the same architecture to other detectors (e.g., **CLD**), or

* you wish to include additional physics processes (e.g., background).

To proceed, you need to clone [this repository](https://github.com/andread3vita/Tracking_DC/tree/devBranch), which provides the Python implementation of the model along with instructions for re-training it and converting the resulting checkpoint file (`.ckpt`) into an **ONNX** file. The ONNX file can then be used to load the model for inference.

## Genfit Track Fitter

**Genfit Track Fitter** is a Gaudi MultiTransformer algorithm developed within the Key4hep framework. It refines reconstructed tracks using the GENFIT tracking toolkit.

The algorithm performs a full track fit including:

- Magnetic field propagation
- Material effects handling
- Drift chamber and planar measurement support

It outputs a new edm4hep::TrackCollection containing fitted tracks with updated track states and fit quality information.

### Technical Implementation

The **GenFit-based Track Fitter** has been integrated into the reconstruction chain by developing a dedicated interface layer between **EDM4hep/Key4hep** data structures and the **GenFit** tracking framework.

The main fitting algorithm is implemented in:

- `components/GenfitTrackFitter.cpp`

while the interface layer connecting EDM4hep and GenFit is located in:

- `include/genfit_interfaces/`
- `src/genfit_interfaces/`

#### 1. EDM4hep - GenFit Interface Layer

To ensure a clean separation between the experiment data model and the fitting engine, dedicated wrapper classes were implemented.

##### Implemented Interfaces

- `GenfitTrack.{h,cpp}`
  Converts an `edm4hep::Track` into a GenFit track object and back.

- `GenfitPlanarMeasurement.{h,cpp}`
  Wraps silicon planar hits into GenFit `PlanarMeasurement` objects.

- `GenfitWireMeasurement.{h,cpp}`
  Converts drift chamber hits into GenFit `WireMeasurement` objects.

- `GenfitField.{h,cpp}`
  Provides the magnetic field interface required by GenFit.

- `GenfitMaterialInterface.{h,cpp}`
  Connects detector material effects (energy loss, multiple scattering) to GenFit.

This modular structure ensures:
- EDM4hep objects remain unchanged.
- GenFit operates with its native abstractions.
- All conversion logic is centralized and reusable.

#### 2. Track Preparation

Inside `GenfitTrackFitter.cpp`, the workflow proceeds as follows:

1. Read the input `edm4hep::TrackCollection`.
2. Extract associated tracker hits.
3. Convert hits into:
   - `GenfitPlanarMeasurement` (silicon detectors),
   - `GenfitWireMeasurement` (drift chamber).
4. Build a `GenfitTrack` object containing all measurements.

#### 3. Magnetic Field and Material Setup

Before running the fit, GenFit requires:

- A magnetic field instance provided by `GenfitField`
- A material effects interface provided by `GenfitMaterialInterface`

These are initialized and registered within GenFit’s global environment.

#### 4. Track Fitting

The fitting procedure consists of:

1. Instantiating a Kalman-based fitter (e.g. `KalmanFitterRefTrack`).
2. Executing the fit on the `GenfitTrack`.
3. Extracting fitted parameters:
   - Position
   - Momentum
   - Covariance matrix
   - $\chi^2$ and fit quality indicators

#### 5. GenFit - EDM4hep Conversion

After the fit:

1. Fitted parameters are converted back to EDM4hep format.
2. A new `edm4hep::Track` is created.
3. Track states and covariance matrices are stored.
4. The results are written to the output collection.


