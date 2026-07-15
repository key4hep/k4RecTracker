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
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// ONNX Runtime and Torch
#include "onnxruntime_cxx_api.h"
#include <torch/torch.h>

// ROOT
#include "TVector3.h"

// Gaudi
#include "Gaudi/Property.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4HEP
#include "edm4hep/SenseWireHitCollection.h"
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
 * @note Cluster identifier zero is treated as noise and does not create a track.
 */
struct GGTFTrackFinder final : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection>(
                                   const std::vector<const edm4hep::TrackerHitPlaneCollection*>&,
                                   const std::vector<const edm4hep::SenseWireHitCollection*>&)> {

  GGTFTrackFinder(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {KeyValues("InputPlanarHitCollections", {"InputPlanarHitCollections"}),
                          KeyValues("InputWireHitCollections", {"InputWireHitCollections"})},
                         {KeyValues("OutputTracksGGTF", {"OutputTracksGGTF"})}) {}

  StatusCode initialize() override {
    if (m_modelPath.value().empty()) {
      error() << "ModelPath is empty; an ONNX model file must be configured." << endmsg;
      return StatusCode::FAILURE;
    }

    try {
      m_memoryInfo =
          std::make_unique<Ort::MemoryInfo>(Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault));
      m_environment = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "GGTFTrackFinder");

      m_sessionOptions.SetIntraOpNumThreads(1);
      m_sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

      m_session = std::make_unique<Ort::Session>(*m_environment, m_modelPath.value().c_str(), m_sessionOptions);

      if (m_session->GetInputCount() != 1 || m_session->GetOutputCount() != 1) {
        error() << "The configured model must expose exactly one input and one output." << endmsg;
        return StatusCode::FAILURE;
      }

      Ort::AllocatorWithDefaultOptions allocator;
      m_inputName = m_session->GetInputNameAllocated(0, allocator).get();
      m_outputName = m_session->GetOutputNameAllocated(0, allocator).get();
    } catch (const Ort::Exception& exception) {
      error() << "Failed to initialize ONNX Runtime: " << exception.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (const std::exception& exception) {
      error() << "Failed to initialize GGTFTrackFinder: " << exception.what() << endmsg;
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::TrackCollection>
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& inputPlanarHitCollections,
             const std::vector<const edm4hep::SenseWireHitCollection*>& inputWireHitCollections) const override {
    const std::uint64_t eventNumber = m_eventCounter.fetch_add(1, std::memory_order_relaxed);
    info() << "Processing event number: " << eventNumber << endmsg;

    edm4hep::TrackCollection outputTracks;

    HitBatch batch;
    reserveBatch(inputPlanarHitCollections, inputWireHitCollections, batch);
    appendPlanarHits(inputPlanarHitCollections, batch);
    appendWireHits(inputWireHitCollections, batch);

    if (batch.nHits == 0) {
      return std::make_tuple(std::move(outputTracks));
    }

    if (batch.nHits > kMaxHits) {
      warning() << "Event " << eventNumber << " has " << batch.nHits << " hits, exceeding the configured limit of "
                << kMaxHits << "; skipping." << endmsg;
      return std::make_tuple(std::move(outputTracks));
    }

    const std::vector<float> modelOutput = runInference(batch.features, batch.nHits);
    const torch::Tensor clusterIds = get_clustering(modelOutput, batch.nHits, m_tbeta, m_td);

    buildTracks(clusterIds, batch, inputPlanarHitCollections, inputWireHitCollections, outputTracks);

    return std::make_tuple(std::move(outputTracks));
  }

  StatusCode finalize() override {
    info() << "Run report:" << endmsg;
    info() << "Number of analysed events: " << m_eventCounter.load(std::memory_order_relaxed) << endmsg;
    info() << "----------------" << endmsg;
    return StatusCode::SUCCESS;
  }

private:
  enum class HitType : std::uint8_t { Planar = 0, Wire = 1 };

  /**
   * @struct HitBatch
   * @brief Flattened per-event feature buffer and source-hit bookkeeping.
   */
  struct HitBatch {
    std::vector<float> features;                // Seven float features per hit.
    std::vector<HitType> hitTypes;              // Original hit kind.
    std::vector<std::size_t> collectionIndices; // Source collection index.
    std::vector<std::size_t> hitIndices;        // Index within the source collection.
    std::size_t nHits = 0;                      // Total number of flattened hits.
  };

  /**
   * @brief Reserve enough storage for all hits in the event.
   *
   * Reserving all parallel arrays once avoids repeated reallocations and keeps
   * their capacities consistent.
   *
   * @param planarCollections Planar hit collections.
   * @param wireCollections Wire hit collections.
   * @param batch Batch whose buffers will be reserved.
   */
  static void reserveBatch(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& planarCollections,
                           const std::vector<const edm4hep::SenseWireHitCollection*>& wireCollections,
                           HitBatch& batch) {
    std::size_t totalHits = 0;
    for (const auto* collection : planarCollections) {
      if (collection != nullptr) {
        totalHits += collection->size();
      }
    }
    for (const auto* collection : wireCollections) {
      if (collection != nullptr) {
        totalHits += collection->size();
      }
    }

    if (totalHits > static_cast<std::size_t>(kMaxHits)) {
      // Reserve only up to the accepted event size. The event will be rejected
      // after flattening reaches the configured limit.
      totalHits = static_cast<std::size_t>(kMaxHits) + 1U;
    }

    batch.features.reserve(totalHits * kFeatureCount);
    batch.hitTypes.reserve(totalHits);
    batch.collectionIndices.reserve(totalHits);
    batch.hitIndices.reserve(totalHits);
  }

  /**
   * @brief Append planar tracker hits to the flattened model input.
   *
   * A planar hit is encoded as `(x, y, z, 1, 0, 0, 0)`.
   *
   * @param collections Input planar hit collections.
   * @param batch Destination event batch.
   */
  static void appendPlanarHits(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& collections,
                               HitBatch& batch) {
    for (std::size_t collectionIndex = 0; collectionIndex < collections.size(); ++collectionIndex) {
      const auto* collection = collections[collectionIndex];
      if (collection == nullptr) {
        throw std::invalid_argument("Null planar-hit collection pointer");
      }

      for (std::size_t hitIndex = 0; hitIndex < collection->size(); ++hitIndex) {
        if (batch.nHits > static_cast<std::size_t>(kMaxHits)) {
          return;
        }

        const auto hit = collection->at(hitIndex);
        const auto position = hit.getPosition();

        batch.features.insert(batch.features.end(), {static_cast<float>(position.x), static_cast<float>(position.y),
                                                     static_cast<float>(position.z), 1.0F, 0.0F, 0.0F, 0.0F});
        batch.hitTypes.push_back(HitType::Planar);
        batch.collectionIndices.push_back(collectionIndex);
        batch.hitIndices.push_back(hitIndex);
        ++batch.nHits;
      }
    }
  }

  /**
   * @brief Append wire hits to the flattened model input.
   *
   * Each drift-chamber hit is represented by one left/right ambiguity segment:
   * the first three features contain the left point, the fourth identifies a
   * wire hit, and the final three contain the vector from left to right.
   *
   * @param collections Input wire-hit collections.
   * @param batch Destination event batch.
   */
  static void appendWireHits(const std::vector<const edm4hep::SenseWireHitCollection*>& collections, HitBatch& batch) {
    for (std::size_t collectionIndex = 0; collectionIndex < collections.size(); ++collectionIndex) {
      const auto* collection = collections[collectionIndex];
      if (collection == nullptr) {
        throw std::invalid_argument("Null wire-hit collection pointer");
      }

      for (std::size_t hitIndex = 0; hitIndex < collection->size(); ++hitIndex) {
        if (batch.nHits > static_cast<std::size_t>(kMaxHits)) {
          return;
        }

        const auto hit = collection->at(hitIndex);
        const edm4hep::Vector3d wirePosition = hit.getPosition();
        const TVector3 wirePositionVector(wirePosition.x, wirePosition.y, wirePosition.z);

        const double distanceToWire = hit.getDistanceToWire();
        const double azimuthalAngle = hit.getWireAzimuthalAngle();
        const double stereoAngle = hit.getWireStereoAngle();

        if (!std::isfinite(distanceToWire) || !std::isfinite(azimuthalAngle) || !std::isfinite(stereoAngle)) {
          throw std::runtime_error("Wire hit contains non-finite geometry values");
        }

        TVector3 zPrime, xPrime, yPrime;
        const double dx = std::sin(stereoAngle) * std::sin(azimuthalAngle);
        const double dy = -std::sin(stereoAngle) * std::cos(azimuthalAngle);
        const double dz = std::cos(stereoAngle);
        zPrime = TVector3(dx, dy, dz).Unit();
        xPrime = TVector3(1.0, 0.0, -dx / dz).Unit(); // x' = normalize([1, 0, -dx/dz])
        yPrime = zPrime.Cross(xPrime).Unit();         // y' = z' x x'

        const TVector3 leftLocal(-distanceToWire, 0.0, 0.0);
        const TVector3 rightLocal(distanceToWire, 0.0, 0.0);

        const TVector3 leftGlobal =
            xPrime * leftLocal.X() + yPrime * leftLocal.Y() + zPrime * leftLocal.Z() + wirePositionVector;
        const TVector3 rightGlobal =
            xPrime * rightLocal.X() + yPrime * rightLocal.Y() + zPrime * rightLocal.Z() + wirePositionVector;
        const TVector3 ambiguityVector = rightGlobal - leftGlobal;

        batch.features.insert(batch.features.end(),
                              {static_cast<float>(leftGlobal.X()), static_cast<float>(leftGlobal.Y()),
                               static_cast<float>(leftGlobal.Z()), 0.0F, static_cast<float>(ambiguityVector.X()),
                               static_cast<float>(ambiguityVector.Y()), static_cast<float>(ambiguityVector.Z())});
        batch.hitTypes.push_back(HitType::Wire);
        batch.collectionIndices.push_back(collectionIndex);
        batch.hitIndices.push_back(hitIndex);
        ++batch.nHits;
      }
    }
  }

  /**
   * @brief Run the ONNX model for one event.
   *
   * @param features Flattened input feature buffer; must contain `nHits * 7`
   *        float values.
   * @param nHits Number of hits represented by the input buffer.
   * @return Flat model output containing four float values per hit.
   * @throws std::logic_error if the session is not initialized.
   * @throws std::runtime_error if the input or output shape is invalid.
   */
  std::vector<float> runInference(std::vector<float>& features, std::size_t nHits) const {

    if (m_session == nullptr || m_memoryInfo == nullptr) {
      throw std::logic_error("ONNX Runtime session is not initialized");
    }
    if (features.size() != nHits * kFeatureCount) {
      throw std::runtime_error("Feature buffer size does not match nHits * 7");
    }
    if (nHits > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
      throw std::overflow_error("Hit count cannot be represented by the ONNX tensor shape");
    }

    const std::vector<std::int64_t> inputShape = {static_cast<std::int64_t>(nHits),
                                                  static_cast<std::int64_t>(kFeatureCount)};

    Ort::Value input = Ort::Value::CreateTensor<float>(*m_memoryInfo, features.data(), features.size(),
                                                       inputShape.data(), inputShape.size());

    const char* inputNames[] = {m_inputName.c_str()};
    const char* outputNames[] = {m_outputName.c_str()};
    auto outputs = m_session->Run(Ort::RunOptions{nullptr}, inputNames, &input, 1, outputNames, 1);

    if (outputs.size() != 1 || !outputs.front().IsTensor()) {
      throw std::runtime_error("ONNX model did not return exactly one tensor output");
    }

    const auto outputInfo = outputs.front().GetTensorTypeAndShapeInfo();
    const std::size_t outputElementCount = outputInfo.GetElementCount();
    const std::size_t expectedElementCount = nHits * kOutputValuesPerHit;
    if (outputElementCount != expectedElementCount) {
      throw std::runtime_error("Unexpected ONNX output size: expected four values per hit");
    }

    const float* outputData = outputs.front().GetTensorData<float>();
    return {outputData, outputData + outputElementCount};
  }

  /**
   * @brief Group hits by cluster identifier and emit one track per cluster.
   *
   * @param clusterIds One-dimensional tensor with one cluster identifier per hit.
   * @param batch Flattened hit metadata used to recover source EDM objects.
   * @param planarCollections Original planar hit collections.
   * @param wireCollections Original wire hit collections.
   * @param outputTracks Destination track collection.
   */
  static void buildTracks(const torch::Tensor& clusterIds, const HitBatch& batch,
                          const std::vector<const edm4hep::TrackerHitPlaneCollection*>& planarCollections,
                          const std::vector<const edm4hep::SenseWireHitCollection*>& wireCollections,
                          edm4hep::TrackCollection& outputTracks) {
    if (!clusterIds.defined()) {
      throw std::runtime_error("Clustering returned an undefined tensor");
    }

    const torch::Tensor ids = clusterIds.to(torch::kCPU).to(torch::kInt64).contiguous().view({-1});
    if (static_cast<std::size_t>(ids.numel()) != batch.nHits) {
      throw std::runtime_error("Cluster-id count does not match the number of input hits");
    }

    const std::int64_t nHits = ids.numel();
    if (nHits == 0) {
      return;
    }

    const torch::Tensor order = torch::argsort(ids);
    const torch::Tensor sortedIds = ids.index_select(0, order);
    const auto orderAccessor = order.accessor<std::int64_t, 1>();
    const auto sortedIdsAccessor = sortedIds.accessor<std::int64_t, 1>();

    std::int64_t clusterBegin = 0;
    while (clusterBegin < nHits) {
      std::int64_t clusterEnd = clusterBegin + 1;
      while (clusterEnd < nHits && sortedIdsAccessor[clusterEnd] == sortedIdsAccessor[clusterBegin]) {
        ++clusterEnd;
      }

      const std::int64_t clusterId = sortedIdsAccessor[clusterBegin];
      if (clusterId > 0) {

        auto outputTrack = outputTracks.create();
        outputTrack.setType(1);

        for (std::int64_t index = clusterBegin; index < clusterEnd; ++index) {
          const std::size_t row = static_cast<std::size_t>(orderAccessor[index]);
          const std::size_t collectionIndex = batch.collectionIndices[row];
          const std::size_t hitIndex = batch.hitIndices[row];

          if (batch.hitTypes[row] == HitType::Planar) {
            outputTrack.addToTrackerHits(planarCollections[collectionIndex]->at(hitIndex));
          } else {
            outputTrack.addToTrackerHits(wireCollections[collectionIndex]->at(hitIndex));
          }
        }
      }

      clusterBegin = clusterEnd;
    }
  }

  static constexpr std::size_t kFeatureCount = 7;
  static constexpr std::size_t kOutputValuesPerHit = 4;
  static constexpr std::size_t kMaxHits = 20000;

  mutable std::atomic<std::uint64_t> m_eventCounter{0}; // Thread-safe processed-event counter.

  std::unique_ptr<Ort::Env> m_environment;       // ONNX Runtime environment.
  std::unique_ptr<Ort::Session> m_session;       // ONNX inference session.
  Ort::SessionOptions m_sessionOptions;          // ONNX session configuration.
  std::unique_ptr<Ort::MemoryInfo> m_memoryInfo; // CPU tensor-memory descriptor.
  std::string m_inputName;                       // Owned ONNX input name.
  std::string m_outputName;                      // Owned ONNX output name.

  Gaudi::Property<std::string> m_modelPath{this, "ModelPath", "", "Path to the ONNX model file"};
  Gaudi::Property<double> m_tbeta{this, "Tbeta", 0.6, "Threshold used to identify cluster core points"};
  Gaudi::Property<double> m_td{this, "Td", 0.3, "Radius used to assign nearby hits to a cluster core"};
};

DECLARE_COMPONENT(GGTFTrackFinder)