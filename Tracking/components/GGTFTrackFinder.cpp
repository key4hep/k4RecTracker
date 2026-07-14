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
#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include "utils.h"

/** @struct GGTFTrackFinder
 *
 * Gaudi MultiTransformer that generates a Track collection by analyzing the digitalized hits through the
 * GGTFTrackFinder. The first step takes the raw hits and it returns a collection of 4-dimensional points inside an
 * embedding space. Each 4-dim point has 3 geometric coordinates and 1 charge, the meaning of which can be described
 * intuitively by a potential, which attracts hits belonging to the same cluster and drives away those that do not.
 * This collection of 4-dim points is analysed by a clustering step, which groups together hits belonging to the
 * same track.
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

    m_fInfo = std::make_unique<Ort::MemoryInfo>(Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault));
    m_fEnv = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");

    m_fSessionOptions.SetIntraOpNumThreads(1);
    m_fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);

    m_fSession = std::make_unique<Ort::Session>(*m_fEnv, m_modelPath.value().c_str(), m_fSessionOptions);

    Ort::AllocatorWithDefaultOptions allocator;
    m_inputNamesOwned.emplace_back(m_fSession->GetInputNameAllocated(0, allocator).get());
    m_outputNamesOwned.emplace_back(m_fSession->GetOutputNameAllocated(0, allocator).get());
    m_fInames.push_back(m_inputNamesOwned.back().c_str());
    m_fOnames.push_back(m_outputNamesOwned.back().c_str());

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::TrackCollection>
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& inputPlanarHitCollections,
             const std::vector<const edm4hep::SenseWireHitCollection*>& inputWireHitCollections) const override {

    info() << "Processing event number: " << m_indexCounter++ << endmsg;

    /////////////////////////////
    ///// Hit Preprocessing /////
    /////////////////////////////

    HitBatch batch;
    appendPlanarHits(inputPlanarHitCollections, batch);
    appendWireHits(inputWireHitCollections, batch);

    edm4hep::TrackCollection outputTracks;

    if (batch.nHits == 0) {
      return std::make_tuple(std::move(outputTracks));
    }
    if (batch.nHits >= kMaxHits) {
      warning() << "Event " << m_indexCounter << " has " << batch.nHits << " hits, exceeding the configured limit of "
                << kMaxHits << " - skipping." << endmsg;
      return std::make_tuple(std::move(outputTracks));
    }

    ////////////////////////////
    ///// Run Track Finder /////
    ////////////////////////////

    const std::vector<float> modelOutput = runInference(batch.features, batch.nHits);
    const torch::Tensor clusterIds = get_clustering(modelOutput, batch.nHits, m_tbeta, m_td);

    //////////////////////////
    ///// Output Results /////
    //////////////////////////

    buildTracks(clusterIds, batch, inputPlanarHitCollections, inputWireHitCollections, outputTracks);

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
  // Flattened per-event hit buffer
  struct HitBatch {
    std::vector<float> features;          // size = nHits * 7
    std::vector<int64_t> hitType;         // 0 = planar, 1 = wire
    std::vector<int64_t> collectionIndex; // which input collection the hit came from
    std::vector<int64_t> subIndex;        // index of the hit within that collection
    int64_t nHits = 0;
  };

  void appendPlanarHits(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& collections,
                        HitBatch& batch) const {
    int64_t collIdx = 0;
    for (const auto* collection : collections) {
      batch.features.reserve(batch.features.size() + collection->size() * 7);
      int64_t subIdx = 0;
      for (const auto& hit : *collection) {
        const auto pos = hit.getPosition();
        batch.features.push_back(static_cast<float>(pos.x));
        batch.features.push_back(static_cast<float>(pos.y));
        batch.features.push_back(static_cast<float>(pos.z));
        batch.features.push_back(1.0f);
        batch.features.push_back(0.0f);
        batch.features.push_back(0.0f);
        batch.features.push_back(0.0f);

        batch.hitType.push_back(0);
        batch.collectionIndex.push_back(collIdx);
        batch.subIndex.push_back(subIdx);

        ++batch.nHits;
        ++subIdx;
      }
      ++collIdx;
    }
  }

  void appendWireHits(const std::vector<const edm4hep::SenseWireHitCollection*>& collections, HitBatch& batch) const {
    int64_t collIdx = 0;
    for (const auto* collection : collections) {
      batch.features.reserve(batch.features.size() + collection->size() * 7);
      int64_t subIdx = 0;
      for (const auto& hit : *collection) {

        const edm4hep::Vector3d wirePos = hit.getPosition();
        const TVector3 wirePosVector(wirePos.x, wirePos.y, wirePos.z);

        const double distanceToWire = hit.getDistanceToWire();
        const double wireAzimuthalAngle = hit.getWireAzimuthalAngle();
        const double wireStereoAngle = hit.getWireStereoAngle();

        TVector3 zPrime, xPrime, yPrime;
        const double dx = std::sin(wireStereoAngle) * std::sin(wireAzimuthalAngle);
        const double dy = -std::sin(wireStereoAngle) * std::cos(wireAzimuthalAngle);
        const double dz = std::cos(wireStereoAngle);
        zPrime = TVector3(dx, dy, dz).Unit();
        xPrime = TVector3(1.0, 0.0, -dx / dz).Unit(); // x' = normalize([1, 0, -dx/dz])
        yPrime = zPrime.Cross(xPrime).Unit();         // y' = z' x x'

        const TVector3 leftLocal(-distanceToWire, 0.0, 0.0);
        const TVector3 rightLocal(distanceToWire, 0.0, 0.0);

        const TVector3 leftGlobal =
            xPrime * leftLocal.X() + yPrime * leftLocal.Y() + zPrime * leftLocal.Z() + wirePosVector;
        const TVector3 rightGlobal =
            xPrime * rightLocal.X() + yPrime * rightLocal.Y() + zPrime * rightLocal.Z() + wirePosVector;

        batch.features.push_back(static_cast<float>(leftGlobal.X()));
        batch.features.push_back(static_cast<float>(leftGlobal.Y()));
        batch.features.push_back(static_cast<float>(leftGlobal.Z()));
        batch.features.push_back(0.0f);
        batch.features.push_back(static_cast<float>(rightGlobal.X() - leftGlobal.X()));
        batch.features.push_back(static_cast<float>(rightGlobal.Y() - leftGlobal.Y()));
        batch.features.push_back(static_cast<float>(rightGlobal.Z() - leftGlobal.Z()));

        batch.hitType.push_back(1);
        batch.collectionIndex.push_back(collIdx);
        batch.subIndex.push_back(subIdx);

        ++batch.nHits;
        ++subIdx;
      }
      ++collIdx;
    }
  }

  std::vector<float> runInference(std::vector<float>& features, int64_t nHits) const {
    const std::vector<int64_t> inputShape = {nHits, 7};

    std::vector<Ort::Value> inputs;
    inputs.emplace_back(Ort::Value::CreateTensor<float>(*m_fInfo, features.data(), features.size(), inputShape.data(),
                                                        inputShape.size()));

    auto outputs = m_fSession->Run(Ort::RunOptions{nullptr}, m_fInames.data(), inputs.data(), m_fInames.size(),
                                   m_fOnames.data(), m_fOnames.size());

    const float* raw = outputs.front().GetTensorMutableData<float>();
    return std::vector<float>(raw, raw + nHits * 4);
  }

  // Groups hits by cluster id and emits one output track per cluster.
  void buildTracks(const torch::Tensor& clusterIds, const HitBatch& batch,
                   const std::vector<const edm4hep::TrackerHitPlaneCollection*>& planarCollections,
                   const std::vector<const edm4hep::SenseWireHitCollection*>& wireCollections,
                   edm4hep::TrackCollection& outputTracks) const {

    const auto ids = clusterIds.to(torch::kInt64).contiguous();
    const int64_t n = ids.size(0);
    if (n == 0) {
      return;
    }

    const auto order = torch::argsort(ids);
    const auto sortedIds = ids.index_select(0, order);

    const auto orderAcc = order.accessor<int64_t, 1>();
    const auto sortedIdsAcc = sortedIds.accessor<int64_t, 1>();

    int64_t start = 0;
    while (start < n) {
      int64_t end = start + 1;
      while (end < n && sortedIdsAcc[end] == sortedIdsAcc[start]) {
        ++end;
      }

      // Noise/unclustered hits (cluster id 0) do not produce a track.
      if (sortedIdsAcc[start] == 0) {
        start = end;
        continue;
      }

      auto outputTrack = outputTracks.create();
      outputTrack.setType(1);

      for (int64_t k = start; k < end; ++k) {
        const int64_t row = orderAcc[k];
        if (batch.hitType[row] == 0) {
          outputTrack.addToTrackerHits(planarCollections[batch.collectionIndex[row]]->at(batch.subIndex[row]));
        } else {
          outputTrack.addToTrackerHits(wireCollections[batch.collectionIndex[row]]->at(batch.subIndex[row]));
        }
      }

      start = end;
    }
  }

private:
  static constexpr int64_t kMaxHits = 20000;

  // Pointer to the ONNX environment. Manages global state such as logging and threading.
  std::unique_ptr<Ort::Env> m_fEnv;

  // Pointer to the ONNX inference session used to execute the model.
  std::unique_ptr<Ort::Session> m_fSession;

  // ONNX session options.
  Ort::SessionOptions m_fSessionOptions;

  // ONNX memory info, describing memory type/device for tensor creation.
  std::unique_ptr<Ort::MemoryInfo> m_fInfo;

  // Owned storage for input/output names
  std::vector<std::string> m_inputNamesOwned;
  std::vector<std::string> m_outputNamesOwned;
  std::vector<const char*> m_fInames;
  std::vector<const char*> m_fOnames;

  // Property to specify the path to the ONNX model file.
  Gaudi::Property<std::string> m_modelPath{this, "ModelPath", "", "ModelPath"};

  // Properties of the clustering step.
  Gaudi::Property<double> m_tbeta{this, "Tbeta", 0.6,
                                  "tbeta: threshold used to identify core points in clusters (tracks)"};
  Gaudi::Property<double> m_td{this, "Td", 0.3,
                               "td: radius around a core point used to assign nearby hits to its cluster"};
};

DECLARE_COMPONENT(GGTFTrackFinder)