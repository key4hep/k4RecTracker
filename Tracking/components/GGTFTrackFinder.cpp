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
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// ONNX Runtime
#include "onnxruntime_cxx_api.h"

// Gaudi
#include "Gaudi/Property.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// EDM4HEP
#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

/** @struct GGTFTrackFinder
 *
 * Gaudi MultiTransformer that generates a Track collection by analyzing the digitalized hits through the
 * GGTFTrackFinder. The first step takes the raw hits and it returns a collection of points inside an embedding
 * space whose dimension is configured through the EmbeddingCoordinates property, plus one trailing charge value.
 * The meaning of the charge can be described intuitively by a potential, which attracts hits belonging to the
 * same cluster and drives away those that do not. This collection of points is analysed by a clustering step,
 * which groups together hits belonging to the same track.
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
    if (m_embeddingCoordinates.value() == 0) {
      error() << "EmbeddingCoordinates must be greater than 0." << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_inferenceThreads.value() < 1) {
      error() << "InferenceThreads must be at least 1." << endmsg;
      return StatusCode::FAILURE;
    }

    try {
      m_memoryInfo =
          std::make_unique<Ort::MemoryInfo>(Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault));
      m_environment = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "GGTFTrackFinder");

      m_sessionOptions.SetIntraOpNumThreads(m_inferenceThreads.value());
      // Lower optimization levels prevent operator fusion and constant folding, which noticeably
      // slows down inference for the matmul-heavy GATr model.
      m_sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
      // The memory arena stays disabled by design: without it the large per-event intermediate
      // tensors are released to the system allocator when the run ends, keeping the held heap at
      // tens of MB instead of GBs. SessionRecreateEvents caps the residual per-shape retention.
      m_sessionOptions.DisableCpuMemArena();
      // The number of hits changes from event to event and the memory-pattern machinery records an
      // allocation plan for every distinct input shape, growing the heap by ~200 KB per event for
      // the lifetime of the job; allocation patterns only pay off for static shapes anyway.
      m_sessionOptions.DisableMemPattern();

      m_session = std::make_shared<Ort::Session>(*m_environment, m_modelPath.value().c_str(), m_sessionOptions);

      if (m_session->GetInputCount() < 1 || m_session->GetOutputCount() < 1) {
        error() << "The configured model must expose at least one input and one output." << endmsg;
        return StatusCode::FAILURE;
      }

      Ort::AllocatorWithDefaultOptions allocator;
      m_inputName = m_session->GetInputNameAllocated(0, allocator).get();
      m_outputName = m_session->GetOutputNameAllocated(0, allocator).get();

      // Models may return extra per-hit values, so only a lower bound on the values per hit is checked.
      const auto outputShape = m_session->GetOutputTypeInfo(0).GetTensorTypeAndShapeInfo().GetShape();
      if (outputShape.size() != 2) {
        error() << "Expected a rank-2 model output [nHits, values per hit], got rank " << outputShape.size() << endmsg;
        return StatusCode::FAILURE;
      }
      // outputShape[1] is -1 if the dimension is dynamic, then it can only be checked for each event
      if (outputShape[1] > 0 && static_cast<std::size_t>(outputShape[1]) < m_embeddingCoordinates.value() + 1) {
        error() << "The model returns " << outputShape[1]
                << " values per hit but at least EmbeddingCoordinates + 1 = " << m_embeddingCoordinates.value() + 1
                << " are needed" << endmsg;
        return StatusCode::FAILURE;
      }
    } catch (const Ort::Exception& exception) {
      error() << "Failed to initialize ONNX Runtime: " << exception.what() << endmsg;
      return StatusCode::FAILURE;
    } catch (const std::exception& exception) {
      error() << "Failed to initialize GGTFTrackFinder: " << exception.what() << endmsg;
      return StatusCode::FAILURE;
    }

    if (m_session == nullptr || m_memoryInfo == nullptr) {
      error() << "ONNX Runtime session is not initialized" << endmsg;
      return StatusCode::FAILURE;
    }

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::TrackCollection>
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& inputPlanarHitCollections,
             const std::vector<const edm4hep::SenseWireHitCollection*>& inputWireHitCollections) const override {
    const std::uint64_t eventNumber = m_eventCounter.fetch_add(1, std::memory_order_relaxed);

    edm4hep::TrackCollection outputTracks;

    const std::size_t inputHitCount = countHits(inputPlanarHitCollections, inputWireHitCollections);

    // MaxHits is used to prevent events with an exceptionally large number of hits
    // from being processed by the model. Very large inputs can cause ONNX Runtime C++ inference
    // to fail or crash because input size directly affects memory usage and computational cost.
    if (inputHitCount == 0 || inputHitCount >= m_maxHits.value()) {
      if (inputHitCount >= m_maxHits.value()) {
        warning() << "Event " << eventNumber << " has " << inputHitCount << " hits, reaching the configured limit of "
                  << m_maxHits.value() << "; skipping." << endmsg;
      }
      return std::make_tuple(std::move(outputTracks));
    }

    // One workspace per processing thread, reused across events: after the largest event the
    // vectors hold enough capacity and these buffers stop allocating on the heap.
    thread_local Workspace workspace;
    workspace.features.resize(inputHitCount * kFeatureCount);
    workspace.hitRefs.clear();
    workspace.hitRefs.reserve(inputHitCount);

    float* featureOut = workspace.features.data();
    featureOut = appendPlanarHits(inputPlanarHitCollections, featureOut, workspace);
    featureOut = appendWireHits(inputWireHitCollections, featureOut, workspace);

    // Wire hits with non-finite geometry are skipped, so the buffer may be shorter than reserved.
    const std::size_t nHits = workspace.hitRefs.size();
    workspace.features.resize(nHits * kFeatureCount);

    // ORT retains per-input-shape allocations for every distinct hit count it meets (a few
    // MB per new shape for this graph); recreating the session on a fixed cadence releases
    // that state and puts a hard bound on heap growth for jobs processing huge event counts.
    const std::size_t recreateEvery = m_sessionRecreateEvents.value();
    if (recreateEvery > 0 && eventNumber > 0 && eventNumber % recreateEvery == 0) {
      std::lock_guard<std::mutex> sessionLock(m_sessionMutex);
      info() << "Recreating the ONNX session at event " << eventNumber << " (every " << recreateEvery << " events)."
             << endmsg;
      m_session = std::make_shared<Ort::Session>(*m_environment, m_modelPath.value().c_str(), m_sessionOptions);
    }
    std::shared_ptr<Ort::Session> session;
    {
      std::lock_guard<std::mutex> sessionLock(m_sessionMutex);
      session = m_session;
    }

    InferenceOutput inference = runInference(*session, workspace.features);
    clusterHits(inference, workspace);

    buildTracks(workspace.clusterIds, workspace, inputPlanarHitCollections, inputWireHitCollections, outputTracks);

    return std::make_tuple(std::move(outputTracks));
  }

  StatusCode finalize() override {
    info() << "Number of analysed events: " << m_eventCounter.load(std::memory_order_relaxed) << endmsg;
    return StatusCode::SUCCESS;
  }

private:
  enum class HitType : std::uint8_t { Planar = 0, Wire = 1 };

  /**
   * @struct HitRef
   * @brief Bookkeeping needed to recover the source EDM hit of one model input row.
   */
  struct HitRef {
    HitType type;
    std::uint32_t collectionIndex; // Source collection index.
    std::uint32_t hitIndex;        // Index within the source collection.
  };

  /**
   * @struct Workspace
   * @brief Per-event feature buffer, source-hit bookkeeping, and clustering scratch space.
   *
   * One instance per processing thread is reused across events (see the thread_local object in the
   * event loop): the vectors keep their capacity, so after the largest event has been processed the
   * event loop performs no further heap allocations in the feature-extraction and clustering stages.
   */
  struct Workspace {
    std::vector<float> features;           // Seven float features per hit.
    std::vector<HitRef> hitRefs;           // Source-hit bookkeeping, one entry per feature row.
    std::vector<std::int32_t> clusterIds;  // One cluster identifier per hit; 0 means noise.
    std::vector<std::uint32_t> candidates; // Cluster-core candidates in decreasing beta order.
    std::vector<std::uint32_t> unassigned; // Hits not yet absorbed by any cluster.
    std::vector<std::uint32_t> remaining;  // Scratch for the hits surviving one assignment sweep.
    std::vector<std::uint8_t> assigned;    // One flag per hit, nonzero once clustered.
    std::vector<std::uint32_t> counts;     // Hits per cluster identifier during track building.
    std::vector<std::uint32_t> offsets;    // Start row per cluster in the ordered hit list.
    std::vector<std::uint32_t> cursor;     // Per-cluster fill cursor for the ordered hit list.
    std::vector<std::uint32_t> order;      // Hit rows grouped by cluster, input order kept.
  };

  /**
   * @struct InferenceOutput
   * @brief Owns the ONNX output tensor so hits can be clustered in place with no copy.
   *
   * The tensor row stride (valuesPerHit) can exceed EmbeddingCoordinates + 1 when the model
   * appends extra per-hit values; only the leading EmbeddingCoordinates + 1 values are used.
   */
  struct InferenceOutput {
    Ort::Value tensor{nullptr};
    std::size_t nHits = 0;
    std::size_t valuesPerHit = 0; // Row stride of the underlying float buffer.
  };

  /**
   * @brief Count all hits in the event.
   *
   * @param planarCollections Planar hit collections.
   * @param wireCollections Wire hit collections.
   * @return Total number of input hits.
   */
  static std::size_t countHits(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& planarCollections,
                               const std::vector<const edm4hep::SenseWireHitCollection*>& wireCollections) {
    std::size_t totalHits = 0;
    for (const auto* collection : planarCollections) {
      totalHits += collection->size();
    }
    for (const auto* collection : wireCollections) {
      totalHits += collection->size();
    }
    return totalHits;
  }

  /**
   * @brief Append planar tracker hits to the flattened model input.
   *
   * A planar hit is encoded as `(x, y, z, 1, 0, 0, 0)`.
   *
   * @param collections Input planar hit collections.
   * @param featureOut Write position inside the preallocated feature buffer.
   * @param workspace Event workspace (only the source-hit bookkeeping).
   * @return Advanced write position, one row past the last appended hit.
   */
  static float* appendPlanarHits(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& collections,
                                 float* featureOut, Workspace& workspace) {
    std::uint32_t collectionIndex = 0;
    for (const auto* collection : collections) {
      for (std::size_t hitIndex = 0; hitIndex < collection->size(); ++hitIndex) {
        const auto hit = collection->at(hitIndex);
        const auto position = hit.getPosition();
        // Mirror the Python data pipeline: hits with non-finite positions are dropped before inference.
        if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z)) {
          continue;
        }

        featureOut[0] = static_cast<float>(position.x);
        featureOut[1] = static_cast<float>(position.y);
        featureOut[2] = static_cast<float>(position.z);
        featureOut[3] = 1.0F;
        featureOut[4] = 0.0F;
        featureOut[5] = 0.0F;
        featureOut[6] = 0.0F;
        featureOut += kFeatureCount;

        workspace.hitRefs.push_back({HitType::Planar, collectionIndex, static_cast<std::uint32_t>(hitIndex)});
      }
      ++collectionIndex;
    }
    return featureOut;
  }

  /**
   * @brief Append wire hits to the flattened model input.
   *
   * Each drift-chamber hit is represented by one left/right ambiguity segment:
   * the first three features contain the left point, the fourth identifies a
   * wire hit, and the final three contain the vector from left to right.
   *
   * With `x' = normalize(1, 0, -dx/dz)` in the wire-local frame, the left point is
   * `position - d * x'` and the left-to-right vector is `2 * d * x'`, where `d` is the
   * distance to the wire. Only x' enters these expressions, so the features are
   * computed with scalar operations only.
   *
   * @param collections Input wire-hit collections.
   * @param featureOut Write position inside the preallocated feature buffer.
   * @param workspace Event workspace (only the source-hit bookkeeping).
   * @return Advanced write position, one row past the last appended hit.
   */
  float* appendWireHits(const std::vector<const edm4hep::SenseWireHitCollection*>& collections, float* featureOut,
                        Workspace& workspace) const {
    std::uint32_t collectionIndex = 0;
    for (const auto* collection : collections) {
      for (std::size_t hitIndex = 0; hitIndex < collection->size(); ++hitIndex) {
        const auto hit = collection->at(hitIndex);
        const edm4hep::Vector3d wirePosition = hit.getPosition();

        const double distanceToWire = hit.getDistanceToWire();
        const double azimuthalAngle = hit.getWireAzimuthalAngle();
        const double stereoAngle = hit.getWireStereoAngle();

        if (!std::isfinite(distanceToWire) || !std::isfinite(azimuthalAngle) || !std::isfinite(stereoAngle)) {
          warning() << "Skipping wire hit with non-finite geometry values" << endmsg;
          continue;
        }

        // x' = normalize(1, 0, -dx/dz) with dx = sin(stereo) * sin(azimuthal), dz = cos(stereo)
        const double dxOverDz = std::sin(stereoAngle) * std::sin(azimuthalAngle) / std::cos(stereoAngle);
        const double xPrimeNorm = 1.0 / std::sqrt(1.0 + dxOverDz * dxOverDz);
        const double xPrimeX = xPrimeNorm;
        const double xPrimeZ = -dxOverDz * xPrimeNorm;

        featureOut[0] = static_cast<float>(wirePosition.x - distanceToWire * xPrimeX);
        featureOut[1] = static_cast<float>(wirePosition.y); // x' has no y component
        featureOut[2] = static_cast<float>(wirePosition.z - distanceToWire * xPrimeZ);
        featureOut[3] = 0.0F;
        featureOut[4] = static_cast<float>(2.0 * distanceToWire * xPrimeX);
        featureOut[5] = 0.0F;
        featureOut[6] = static_cast<float>(2.0 * distanceToWire * xPrimeZ);
        featureOut += kFeatureCount;

        workspace.hitRefs.push_back({HitType::Wire, collectionIndex, static_cast<std::uint32_t>(hitIndex)});
      }
      ++collectionIndex;
    }
    return featureOut;
  }

  /**
   * @brief Run the ONNX model for one event.
   *
   * The feature buffer is wrapped without copying and the returned object keeps the
   * model output alive so it can be clustered in place.
   *
   * @param features Flattened input feature buffer containing seven float values per hit.
   * @return Inference output owning the model result tensor; nHits is zero for an empty input.
   * @throws std::logic_error if the session is not initialized.
   * @throws std::runtime_error if the input or output shape is invalid.
   */
  InferenceOutput runInference(Ort::Session& session, const std::vector<float>& features) const {
    if (features.size() % kFeatureCount != 0) {
      throw std::runtime_error("Feature buffer size is not a multiple of the feature count");
    }
    const std::size_t nHits = features.size() / kFeatureCount;
    if (nHits > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
      throw std::overflow_error("Hit count cannot be represented by the ONNX tensor shape");
    }
    if (nHits == 0) {
      return InferenceOutput{};
    }

    const std::vector<std::int64_t> inputShape = {static_cast<std::int64_t>(nHits),
                                                  static_cast<std::int64_t>(kFeatureCount)};

    Ort::Value input = Ort::Value::CreateTensor<float>(*m_memoryInfo, const_cast<float*>(features.data()),
                                                       features.size(), inputShape.data(), inputShape.size());

    const char* inputNames[] = {m_inputName.c_str()};
    const char* outputNames[] = {m_outputName.c_str()};
    auto outputs = session.Run(Ort::RunOptions{nullptr}, inputNames, &input, 1, outputNames, 1);

    if (outputs.size() != 1 || !outputs.front().IsTensor()) {
      throw std::runtime_error("ONNX model did not return exactly one tensor output");
    }

    const auto outputInfo = outputs.front().GetTensorTypeAndShapeInfo();
    const std::vector<std::int64_t> outputDims = outputInfo.GetShape();
    const std::size_t valuesPerHit = m_embeddingCoordinates.value() + 1; // embedding coordinates + cluster potential
    // Models may append extra per-hit values after the cluster potential; only the
    // leading embedding coordinates and the cluster potential are used.
    if (outputDims.size() != 2 || outputDims[0] != static_cast<std::int64_t>(nHits)) {
      throw std::runtime_error("Unexpected ONNX output shape: expected [" + std::to_string(nHits) +
                               ", values per hit]");
    }
    const std::size_t outputValuesPerHit = static_cast<std::size_t>(outputDims[1]);
    if (outputValuesPerHit < valuesPerHit) {
      throw std::runtime_error("Unexpected ONNX output size: expected at least " + std::to_string(valuesPerHit) +
                               " values per hit (" + std::to_string(m_embeddingCoordinates.value()) +
                               " embedding coordinates + 1 cluster potential), got " +
                               std::to_string(outputValuesPerHit));
    }

    return InferenceOutput{std::move(outputs.front()), nHits, outputValuesPerHit};
  }

  /**
   * @brief Assign each hit to a cluster based on the model output.
   *
   * Object-condensation style clustering: the model emits a beta logit per hit (the value right
   * after the embedding coordinates); the sigmoid of the logit is the cluster potential, and hits
   * whose potential exceeds Tbeta are cluster core candidates. In decreasing beta
   * order, each unassigned candidate attracts all still-unassigned hits within distance Td in
   * the embedding space. Hits never absorbed keep cluster identifier 0 (noise).
   *
   * @param inference Inference output holding the model result tensor alive.
   * @param workspace Event workspace; its clusterIds member is filled with one cluster
   *                  identifier per hit and cluster identifiers start at 1.
   */
  void clusterHits(const InferenceOutput& inference, Workspace& workspace) const {
    workspace.clusterIds.assign(inference.nHits, std::int32_t{0});
    if (inference.nHits == 0) {
      return;
    }

    const float* output = inference.tensor.GetTensorData<float>();
    const std::size_t stride = inference.valuesPerHit;
    const std::size_t embDims = m_embeddingCoordinates.value();
    const float tbeta = static_cast<float>(m_tbeta.value());
    const float td2 = static_cast<float>(m_td.value()) * static_cast<float>(m_td.value());

    const auto betaAt = [output, stride, embDims](std::uint32_t hitIndex) -> float {
      const float logit = output[static_cast<std::size_t>(hitIndex) * stride + embDims];
      return 1.0F / (1.0F + std::exp(-logit));
    };

    std::vector<std::uint32_t>& candidates = workspace.candidates;
    candidates.clear();
    for (std::uint32_t hitIndex = 0; hitIndex < inference.nHits; ++hitIndex) {
      if (betaAt(hitIndex) > tbeta) {
        candidates.push_back(hitIndex);
      }
    }
    std::stable_sort(candidates.begin(), candidates.end(),
                     [&betaAt](std::uint32_t a, std::uint32_t b) { return betaAt(a) > betaAt(b); });

    workspace.assigned.assign(inference.nHits, std::uint8_t{0});
    std::vector<std::uint8_t>& assigned = workspace.assigned;
    std::vector<std::uint32_t>& unassigned = workspace.unassigned;
    unassigned.resize(inference.nHits);
    std::iota(unassigned.begin(), unassigned.end(), std::uint32_t{0});
    std::vector<std::uint32_t>& remaining = workspace.remaining;
    remaining.clear();

    std::int32_t clusterId = 0;
    std::vector<std::int32_t>& clusterIds = workspace.clusterIds;
    for (const std::uint32_t candidate : candidates) {
      if (unassigned.empty()) {
        break;
      }
      if (assigned[candidate] != 0) {
        continue;
      }

      ++clusterId;
      const float* center = output + static_cast<std::size_t>(candidate) * stride;
      remaining.clear();
      bool anyAssigned = false;
      for (const std::uint32_t hitIndex : unassigned) {
        const float* point = output + static_cast<std::size_t>(hitIndex) * stride;
        float distanceSquared = 0.0F;
        for (std::size_t dim = 0; dim < embDims; ++dim) {
          const float difference = point[dim] - center[dim];
          distanceSquared += difference * difference;
        }
        if (distanceSquared < td2) {
          clusterIds[hitIndex] = clusterId;
          assigned[hitIndex] = 1;
          anyAssigned = true;
        } else {
          remaining.push_back(hitIndex);
        }
      }
      // Mirror the reference implementation: a candidate that attracts nothing (only possible
      // for Td <= 0) is dropped without consuming a cluster identifier.
      if (!anyAssigned) {
        --clusterId;
        continue;
      }
      unassigned.swap(remaining);
    }
  }

  /**
   * @brief Group hits by cluster identifier and emit one track per cluster.
   *
   * Tracks are created in ascending cluster-identifier order and hits within a track keep
   * their original input order. The grouping is a counting sort over the contiguous cluster
   * identifiers produced by clusterHits, reusing the workspace buffers, so no memory is
   * allocated here once the workspace has seen the largest event.
   *
   * @param clusterIds One cluster identifier per hit.
   * @param workspace Event workspace with the source-hit bookkeeping and grouping scratch space.
   * @param planarCollections Original planar hit collections.
   * @param wireCollections Original wire hit collections.
   * @param outputTracks Destination track collection.
   */
  static void buildTracks(const std::vector<std::int32_t>& clusterIds, Workspace& workspace,
                          const std::vector<const edm4hep::TrackerHitPlaneCollection*>& planarCollections,
                          const std::vector<const edm4hep::SenseWireHitCollection*>& wireCollections,
                          edm4hep::TrackCollection& outputTracks) {
    if (clusterIds.size() != workspace.hitRefs.size()) {
      throw std::runtime_error("Cluster-id count does not match the number of input hits");
    }
    if (clusterIds.empty()) {
      return;
    }

    // Cluster ID 0 (and any negative ID) represents noise/unclustered hits and does not produce an output track.
    // Cluster identifiers are contiguous starting at 1, so hits can be grouped with a counting sort.
    std::int32_t maxClusterId = 0;
    for (const std::int32_t clusterId : clusterIds) {
      maxClusterId = std::max(maxClusterId, clusterId);
    }
    if (maxClusterId <= 0) {
      return;
    }

    std::vector<std::uint32_t>& counts = workspace.counts;
    counts.assign(static_cast<std::size_t>(maxClusterId) + 1, std::uint32_t{0});
    for (const std::int32_t clusterId : clusterIds) {
      if (clusterId > 0) {
        ++counts[static_cast<std::size_t>(clusterId)];
      }
    }

    std::vector<std::uint32_t>& offsets = workspace.offsets;
    offsets.resize(counts.size());
    std::uint32_t groupedCount = 0;
    for (std::size_t clusterId = 1; clusterId < counts.size(); ++clusterId) {
      offsets[clusterId] = groupedCount;
      groupedCount += counts[clusterId];
    }

    // Scatter hit rows into one contiguous list grouped by cluster; scanning the rows in
    // input order keeps the hits of every track in their original order.
    std::vector<std::uint32_t>& cursor = workspace.cursor;
    cursor.assign(offsets.begin(), offsets.end());
    std::vector<std::uint32_t>& order = workspace.order;
    order.resize(groupedCount);
    for (std::uint32_t row = 0; row < clusterIds.size(); ++row) {
      if (clusterIds[row] > 0) {
        order[cursor[static_cast<std::size_t>(clusterIds[row])]++] = row;
      }
    }

    for (std::int32_t clusterId = 1; clusterId <= maxClusterId; ++clusterId) {
      if (counts[static_cast<std::size_t>(clusterId)] == 0) {
        continue;
      }
      auto outputTrack = outputTracks.create();
      const std::uint32_t begin = offsets[static_cast<std::size_t>(clusterId)];
      const std::uint32_t end = begin + counts[static_cast<std::size_t>(clusterId)];
      for (std::uint32_t position = begin; position < end; ++position) {
        const HitRef& hitRef = workspace.hitRefs[order[position]];
        if (hitRef.type == HitType::Planar) {
          outputTrack.addToTrackerHits(planarCollections[hitRef.collectionIndex]->at(hitRef.hitIndex));
        } else {
          outputTrack.addToTrackerHits(wireCollections[hitRef.collectionIndex]->at(hitRef.hitIndex));
        }
      }
    }
  }

  static constexpr std::size_t kFeatureCount = 7;
  mutable std::atomic<std::uint64_t> m_eventCounter{0}; // Thread-safe processed-event counter.

  std::unique_ptr<Ort::Env> m_environment;         // ONNX Runtime environment.
  mutable std::shared_ptr<Ort::Session> m_session; // ONNX inference session (swapped on recreation).
  mutable std::mutex m_sessionMutex;               // Serializes session swaps and snapshots.
  Ort::SessionOptions m_sessionOptions;            // ONNX session configuration.
  std::unique_ptr<Ort::MemoryInfo> m_memoryInfo;   // CPU tensor-memory descriptor.
  std::string m_inputName;                         // Owned ONNX input name.
  std::string m_outputName;                        // Owned ONNX output name.

  Gaudi::Property<std::string> m_modelPath{this, "ModelPath", "", "Path to the ONNX model file"};
  Gaudi::Property<std::size_t> m_embeddingCoordinates{
      this, "EmbeddingCoordinates", 3,
      "Number of embedding-space coordinates per hit in the ONNX model output; the model must produce one "
      "additional trailing value per hit, the beta logit (the sigmoid of the logit is the cluster potential "
      "compared against Tbeta)"};
  Gaudi::Property<int> m_inferenceThreads{
      this, "InferenceThreads", 4,
      "Number of intra-operation threads used by the ONNX Runtime session; inference time plateaus beyond "
      "about 4-8 threads for the current GATr model, so increase further only if the job has many free cores"};
  Gaudi::Property<std::size_t> m_maxHits{this, "MaxHits", 20000, "Maximum number of hits accepted per event"};
  Gaudi::Property<std::size_t> m_sessionRecreateEvents{
      this, "SessionRecreateEvents", 0,
      "Recreate the ONNX session every N events (0 disables); bounds the per-input-shape "
      "allocation state ONNX Runtime retains for ever-changing hit counts, at the price of one "
      "session reload per interval"};
  Gaudi::Property<double> m_tbeta{this, "Tbeta", 0.6,
                                  "Threshold on sigmoid(beta logit) used to identify cluster core points"};
  Gaudi::Property<double> m_td{this, "Td", 0.3, "Radius used to assign nearby hits to a cluster core"};
};

DECLARE_COMPONENT(GGTFTrackFinder)
