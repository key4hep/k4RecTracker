/*
 * Copyright (c) 2014-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// Standard Library
#include <algorithm>
#include <cmath>
#include <cstdlib>    
#include <filesystem>
#include <fstream>    
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

// ONNX & Torch
#include "onnxruntime_cxx_api.h"
#include "onnxruntime_run_options_config_keys.h"
#include <ATen/ATen.h>
#include <torch/torch.h>

// ROOT
#include "TVector3.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

// k4FWCore / k4Interface
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// EDM4HEP
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "utils.h"

/** @struct GGTFTrackFinder
 *
 * Gaudi MultiTransformer that generates a Track collection by analyzing the digitalized hits through the
 * GGTFTrackFinder. The first step takes the raw hits and it returns a collection of 4-dimensional points inside an
 * embedding space. Eeach 4-dim point has 3 geometric coordinates and 1 charge, the meaning of which can be described
 * intuitively by a potential, which attracts hits belonging to the same cluster and drives away those that do not. This
 * collection of 4-dim points is analysed by a clustering step, which groups together hits belonging to the same track.
 *
 *  input:
 *    - digitalized hits from DC (global coordinates) : edm4hep::SenseWireHitCollection
 *    - digitalized hits from vertex (global coordinates) : edm4hep::TrackerHitPlaneCollection
 *    - digitalized hits from silicon wrapper (global coordinates) : edm4hep::TrackerHitPlaneCollection
 *
 *  output:
 *    - Track collection : edm4hep::TrackCollection
 *
 *
 *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
 *  @date   2025-11
 *
 */

struct GGTFTrackFinder final : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection>(
                                 const std::vector<const edm4hep::TrackerHitPlaneCollection*>&,
                                 const std::vector<const edm4hep::SenseWireHitCollection*>&)> {

  GGTFTrackFinder(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {

                            KeyValues("InputPlanarHitCollections", {"InputPlanarHitCollections"}),
                            KeyValues("InputWireHitCollections", {"InputWireHitCollections"})

                         },
                         {

                            KeyValues("OutputTracksGGTF", {"OutputTracksGGTF"})

                         }) {}

  StatusCode initialize() override {

    ///////////////////////////////
    ///// ONNX Initialization /////
    ///////////////////////////////

    // Initialize the ONNX memory info object for CPU memory allocation.
    m_fInfo = std::make_unique<Ort::MemoryInfo>(Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault));

    // Create and initialize the ONNX environment with a logging level set to WARNING.
    auto envLocal = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
    m_fEnv = std::move(envLocal);

    // Set the session options to configure the ONNX inference session.
    m_fSessionOptions.SetIntraOpNumThreads(1);

    // Disable all graph optimizations to keep the model execution as close to the original as possible.
    m_fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);

    // Create an ONNX inference session to load the model specified by the `modelPath`.
    auto sessionLocal = std::make_unique<Ort::Session>(*m_fEnv, m_modelPath.value().c_str(), m_fSessionOptions);
    m_fSession = std::move(sessionLocal);

    // Create an ONNX allocator with default options to manage memory allocations during runtime.
    Ort::AllocatorWithDefaultOptions allocator;

    // Get the name of the first input node (index i) from the ONNX model and store it.
    // Get the name of the first output node (index i) from the ONNX model and store it.
    std::size_t i = 0;
    const auto inputNames = m_fSession->GetInputNameAllocated(i, allocator).release();
    const auto outputNames = m_fSession->GetOutputNameAllocated(i, allocator).release();

    m_fInames.push_back(inputNames);
    m_fOnames.push_back(outputNames);

    return StatusCode::SUCCESS;

  }

  std::tuple<edm4hep::TrackCollection>
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& inputPlanarHitCollections,
             const std::vector<const edm4hep::SenseWireHitCollection*>& inputWireHitCollections) const override {

    ////////////////////////////////////////
    ////////// DATA PREPROCESSING //////////
    ////////////////////////////////////////

    info() << "Processing event number: " << m_indexCounter++ << endmsg;

    // Vector to store the global input values for all hits.
    // This will contain position and other hit-specific data to be used as input for the model.
    std::vector<float> listGlobalInputs;
    int globalHitIndex = 0;

    /// Processing hits from the Silicon Detectors
    std::vector<int64_t> listHitTypePlanar;
    int planarHitIndex = 0;
    std::vector<int64_t> listPlanarHitIndices;
    int planarHitCollectionIndex = 0;
    for (const auto& inputHitCollection : inputPlanarHitCollections) {

      int planarHitSubCollectionIndex = 0;
      for (const auto& inputHit : *inputHitCollection) {

        // Add the 3D position of the hit to the global input list.
        listGlobalInputs.push_back(inputHit.getPosition().x);
        listGlobalInputs.push_back(inputHit.getPosition().y);
        listGlobalInputs.push_back(inputHit.getPosition().z);

        // Add placeholder values for additional input dimensions.
        listGlobalInputs.push_back(1.0);
        listGlobalInputs.push_back(0.0);
        listGlobalInputs.push_back(0.0);
        listGlobalInputs.push_back(0.0);

        listHitTypePlanar.push_back(globalHitIndex);

        listPlanarHitIndices.push_back(planarHitCollectionIndex);
        listPlanarHitIndices.push_back(planarHitSubCollectionIndex);

        globalHitIndex += 1;
        planarHitIndex += 1;
        planarHitSubCollectionIndex += 1;
      }

      planarHitCollectionIndex += 1;
    }
    // Convert listHitTypePlanar to a Torch tensor for use in PyTorch models.
    torch::Tensor listHitTypePlanarTensor =
        torch::from_blob(listHitTypePlanar.data(), {planarHitIndex}, torch::kFloat32);

    torch::Tensor listPlanarHitIndices_tensor =
        torch::from_blob(listPlanarHitIndices.data(), {planarHitIndex, 2}, torch::kInt64);


    /// Processing hits from gaseous detectors
    std::vector<int64_t> listHitTypeWire;
    int wireHitIndex = 0;
    std::vector<int64_t> listWireHitIndices;
    int wireHitCollectionIndex = 0;
    for (const auto& inputHitCollection : inputWireHitCollections) {
      int wireHitSubCollectionIndex = 0;
      for (const auto& inputHit : *inputHitCollection) {

        // position along the wire, wire direction and drift distance
        edm4hep::Vector3d wirePos = inputHit.getPosition();
        TVector3 wirePosVector(static_cast<float>(wirePos.x), static_cast<float>(wirePos.y),
                               static_cast<float>(wirePos.z));

        double distanceToWire = inputHit.getDistanceToWire();
        double wireAzimuthalAngle = inputHit.getWireAzimuthalAngle();
        double wireStereoAngle = inputHit.getWireStereoAngle();

        // Direction of the wire: z'
        TVector3 direction(0, 0, 1);
        direction.RotateX(wireStereoAngle);
        direction.RotateZ(wireAzimuthalAngle);

        TVector3 zPrime;
        zPrime = direction.Unit();

        // x' axis, orthogonal to z'
        TVector3 xPrime(1.0, 0.0, -direction.X() / direction.Z());
        xPrime = xPrime.Unit();

        // y' axis = z' × x'
        TVector3 yPrime = zPrime.Cross(xPrime);
        yPrime = yPrime.Unit();

        // Define the local left/right positions
        TVector3 leftLocalPos(-distanceToWire, 0.0, 0.0);
        TVector3 rightLocalPos(distanceToWire, 0.0, 0.0);

        // Transform to global
        TVector3 leftGlobalPos =
            xPrime * leftLocalPos.X() + yPrime * leftLocalPos.Y() + zPrime * leftLocalPos.Z() + wirePosVector;
        TVector3 rightGlobalPos =
            xPrime * rightLocalPos.X() + yPrime * rightLocalPos.Y() + zPrime * rightLocalPos.Z() + wirePosVector;

        // Add the 3D position of the left hit to the global input list.
        listGlobalInputs.push_back(leftGlobalPos.X());
        listGlobalInputs.push_back(leftGlobalPos.Y());
        listGlobalInputs.push_back(leftGlobalPos.Z());

        // Add the difference between the right and left hit positions to the global input list.
        listGlobalInputs.push_back(0.0);
        listGlobalInputs.push_back(rightGlobalPos.X() - leftGlobalPos.X());
        listGlobalInputs.push_back(rightGlobalPos.Y() - leftGlobalPos.Y());
        listGlobalInputs.push_back(rightGlobalPos.Z() - leftGlobalPos.Z());

        listHitTypeWire.push_back(globalHitIndex);

        listWireHitIndices.push_back(wireHitCollectionIndex);
        listWireHitIndices.push_back(wireHitSubCollectionIndex);

        globalHitIndex += 1;
        wireHitIndex += 1;
        wireHitSubCollectionIndex += 1;
      }
    }
    // Convert ListHitType_CDC to a Torch tensor for use in PyTorch models.
    torch::Tensor listHitTypeWireTensor = torch::from_blob(listHitTypeWire.data(), {wireHitIndex}, torch::kFloat32);

    torch::Tensor listWireHitIndicesTensor =
        torch::from_blob(listWireHitIndices.data(), {wireHitIndex, 2}, torch::kInt64);

    torch::Tensor planarTypeTensor = torch::zeros({planarHitIndex, 1}, torch::kInt64);
    torch::Tensor wireTypeTensor = torch::ones({wireHitIndex, 1}, torch::kInt64);

    torch::Tensor planarTypeIndexTensor = torch::cat({planarTypeTensor, listPlanarHitIndices_tensor}, 1);
    torch::Tensor wireTypeIndexTensor = torch::cat({wireTypeTensor, listWireHitIndicesTensor}, 1);
    torch::Tensor listHitIndicesGlobal = torch::cat({planarTypeIndexTensor, wireTypeIndexTensor}, 0);

    //////////////////////////////////
    ////////// TRACK FINDER //////////
    //////////////////////////////////

    // Create a new TrackCollection and TrackerHit3DCollection for storing the output tracks and hits
    edm4hep::TrackCollection outputTracks;
    constexpr int kMaxHits = 20000;
    if (globalHitIndex > 0 && globalHitIndex < kMaxHits) {

      ///////////////////////////////
      ////////// GATR STEP //////////
      ///////////////////////////////

      // Calculate the total size of the input tensor, based on the number of hits (it) and the
      // number of features per hit (7: x, y, z, and four placeholders).
      size_t inputTensorTotalSize = globalHitIndex * 7;
      std::vector<int64_t> inputTensorShape = {globalHitIndex, 7};

      // Create a vector to store the input tensors that will be fed into the ONNX model.
      std::vector<Ort::Value> inputModelTensors;
      inputModelTensors.emplace_back(Ort::Value::CreateTensor<float>(
          *m_fInfo, listGlobalInputs.data(), inputTensorTotalSize, inputTensorShape.data(), inputTensorShape.size()));

      // Run the ONNX inference session with the provided input tensor.
      auto outputModelTensors = m_fSession->Run(Ort::RunOptions{nullptr}, m_fInames.data(), inputModelTensors.data(),
                                                m_fInames.size(), m_fOnames.data(), m_fOnames.size());
      float* floatarr = outputModelTensors.front().GetTensorMutableData<float>();
      std::vector<float> outputModelVector(floatarr, floatarr + globalHitIndex * 4);

      /////////////////////////////////////
      ////////// CLUSTERING STEP //////////
      /////////////////////////////////////

      auto clusteringIndeces = get_clustering(outputModelVector, globalHitIndex, m_tbeta, m_td);
      torch::Tensor uniqueTensor;
      torch::Tensor inverseIndices;
      std::tie(uniqueTensor, inverseIndices) = at::_unique(clusteringIndeces, true, true);

      /////////////////////////////////
      ////////// OUTPUT STEP //////////
      /////////////////////////////////

      // NB: In this implementation, for each track the flag "type" is defined as follows:
      //     0 → unclustered
      //     1 → clustered
      //
      // Here, "unclustered" refers to hits that are not associated with any reconstructed track.
      // This may occur either because the hit originates from background activity or due to
      // an inefficiency/mistake in the track-finding algorithm.

      // Get the total number of unique tracks based on the uniqueTensor size
      int64_t numTracks = uniqueTensor.numel();
      bool has_zero = (uniqueTensor == 0).any().item<bool>();
      if (!has_zero) {
        auto outputTrack = outputTracks.create();
        outputTrack.setType(0);
      }

      // Loop through each unique track ID
      for (int i = 0; i < numTracks; ++i) {

        // Create a new track in the output collection and set its type to the current track ID
        auto outputTrack = outputTracks.create();
        outputTrack.setType(1);

        // Select all hits belonging to the current track
        auto idTrack = uniqueTensor.index({i});
        torch::Tensor mask = (clusteringIndeces == idTrack);
        torch::Tensor indices = torch::nonzero(mask).flatten();

        auto listHitIndicesGlobalView = listHitIndicesGlobal.accessor<int64_t, 2>();
        auto indicesView = indices.accessor<int64_t, 1>();

        for (int64_t j = 0; j < indices.size(0); ++j) {
          int64_t row = indicesView[j];
          int64_t type = listHitIndicesGlobalView[row][0];
          int64_t idx1 = listHitIndicesGlobalView[row][1];
          int64_t idx2 = listHitIndicesGlobalView[row][2];

          if (type == 0) {
            // planar hit
            auto planarCollection = inputPlanarHitCollections[idx1];
            auto hit = planarCollection->at(idx2);
            outputTrack.addToTrackerHits(hit);

          } else {
            // wire hit
            auto wireCollection = inputWireHitCollections[idx1];
            auto hit = wireCollection->at(idx2);
            outputTrack.addToTrackerHits(hit);
          }
        }
      }

      inverseIndices.reset();
      uniqueTensor.reset();
      clusteringIndeces.reset();

      inputModelTensors.clear();
      outputModelTensors.clear();
    }

    listHitTypePlanarTensor.reset();
    listPlanarHitIndices_tensor.reset();
    listHitTypeWireTensor.reset();
    listWireHitIndicesTensor.reset();

    std::vector<int64_t>().swap(listHitTypePlanar);
    std::vector<int64_t>().swap(listPlanarHitIndices);
    std::vector<int64_t>().swap(listHitTypeWire);
    std::vector<int64_t>().swap(listWireHitIndices);

    // Return the output collections as a tuple
    return std::make_tuple(std::move(outputTracks));
  }

  StatusCode finalize() override {

    info() << "Run report:" << endmsg;
    info() << "Number of analysed events: " << m_indexCounter << endmsg;
    info() << "----------------\n" << endmsg;

    return StatusCode::SUCCESS;
  }

public:
  mutable int m_indexCounter = 0;

private:

  // Pointer to the ONNX environment.
  // This object manages the global state of the ONNX runtime, such as logging and threading.
  std::unique_ptr<Ort::Env> m_fEnv;

  // Pointer to the ONNX inference session.
  // This session is used to execute the model for inference.
  std::unique_ptr<Ort::Session> m_fSession;

  // ONNX session options.
  Ort::SessionOptions m_fSessionOptions;

  // ONNX memory info.
  // This object provides information about memory allocation and is used during the creation of
  // ONNX tensors. It specifies the memory type and device (e.g., CPU, GPU).
  std::unique_ptr<Ort::MemoryInfo> m_fInfo;

  // Stores the input and output names for the ONNX model.
  // These vectors contain the names of the inputs (fInames) and outputs (fOnames) that the model expects.
  std::vector<const char*> m_fInames;
  std::vector<const char*> m_fOnames;

  // Property to specify the path to the ONNX model file.
  // This is a configurable property that defines the location of the ONNX model file on the filesystem.
  Gaudi::Property<std::string> m_modelPath{this, "ModelPath", "", "ModelPath"};

  // Properties of the clustering step.
  Gaudi::Property<double> m_tbeta{ this, "Tbeta", 0.6, "tbeta: threshold used to identify core points in clusters (tracks)"};
  Gaudi::Property<double> m_td{this, "Td", 0.3, "td: radius around a core point used to assign nearby hits to its cluster"};


};

DECLARE_COMPONENT(GGTFTrackFinder)