# Key4Hep Tracking Algorithms for FCC-ee detectors

This subfolder contains the implementation of several Tracking tools for FCC-ee detectors using the key4hep framework:

* GGTF_tracking
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

We implemented a **Geometric Graph Track Finder (GGTF)** in `k4RecTracker/Tracking/components/GGTF_tracking.cpp`.  
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
