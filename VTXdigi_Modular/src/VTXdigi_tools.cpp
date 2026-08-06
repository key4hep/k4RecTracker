// VTXdigi_Modular/src/VTXdigi_tools.cpp
#include "VTXdigi_tools.h"
#include <DD4hep/Objects.h>
#include <DD4hep/VolumeManager.h>
#include <edm4hep/Vector3d.h>

namespace VTXdigi_tools {

SimHitWrapper::SimHitWrapper(
  edm4hep::SimTrackerHit simTrackerHit, dd4hep::DDSegmentation::VolumeID volumeID,
  const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder,
  const dd4hep::VolumeManager& volumeManager,
  const std::unique_ptr<dd4hep::rec::CellIDPositionConverter>& cellIDPositionConverter)
    : m_simTrackerHit(simTrackerHit), m_volumeID(volumeID) {

  m_charge = static_cast<float>(m_simTrackerHit.getEDep() * (dd4hep::GeV / dd4hep::keV) * kChargePerkeV); // convert energy deposit (in keV) to number of electrons
  m_layerNumber = GetLayer(m_volumeID, cellIdDecoder);

  // check if the simHit was caused by a primary, secondary or delta particle
  if ( m_simTrackerHit.isProducedBySecondary() ) {
    // if ddsim dropped the MCParticle that caused this simHit, we assume it was a delta ray
    // ddsim drops MCParticles below a certain energy cut to save computing cost and disk space.
    m_mcParticleLevel = MCParticleLevel::Delta;
  }
  else {
    // check if the MCPArticle was created by the generator
    const int32_t simulatorStatus = m_simTrackerHit.getParticle().getSimulatorStatus();
    const int32_t mask = 1 << edm4hep::MCParticle::BITCreatedInSimulation; // should be bit 30
    const bool causedByPrimary = (simulatorStatus & mask) == 0; // bit is not set -> created in generator
    if ( causedByPrimary )
      m_mcParticleLevel = MCParticleLevel::Primary;
    else
    {
      // now check if the MCParticle prod. vertex lies outside this sensors volume (by comparing volumeIDs)
      const edm4hep::Vector3d prodVertex_temp = m_simTrackerHit.getParticle().getVertex();
      const dd4hep::Position prodVertex = 0.1 * ConvertVector_toPosition(prodVertex_temp); // convert edm4hep's mm -> dd4hep's cm
      const dd4hep::DDSegmentation::CellID prodVertex_cellID = cellIDPositionConverter->cellID(prodVertex); // returns 0 if the position is outside of any sensitive volume

      // convert cellID to volumeID - see comment in VTXdigi_Modular::GetVolumeID())
      dd4hep::DDSegmentation::CellID prodVertex_volumeID;
      if (prodVertex_cellID == 0)
        prodVertex_volumeID = 0; // lookupContext(cellID=0) crashes (because cellID 0 does not exist)
      else
        prodVertex_volumeID = volumeManager.lookupContext(prodVertex_cellID)->element.volumeID();

      if (prodVertex_volumeID != m_volumeID) {
        // the MCParticle was created outside of this sensor's sensitive volume
        m_mcParticleLevel = MCParticleLevel::Secondary;
      }
      else {
        m_mcParticleLevel = MCParticleLevel::Delta;
      }
    }
  }
}

void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept {
  std::swap(a.m_simTrackerHit, b.m_simTrackerHit);
  std::swap(a.m_volumeID, b.m_volumeID);
  std::swap(a.m_charge, b.m_charge);
  std::swap(a.m_layerNumber, b.m_layerNumber);
  std::swap(a.m_truthPos, b.m_truthPos);
  std::swap(a.m_mcParticleLevel, b.m_mcParticleLevel);
} // swap(Hit&, Hit&)

// SimulatorStatus bits (see https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_m_c_particle.html)
// 29 : "Backscatter",
// 30 : "CreatedInSimulation",
// 26 : "DecayedInCalorimeter",
// 27 : "DecayedInTracker",
// 22 : "HandledInFastSim",
// 25 : "LeftWorld",
// 23 : "Overlay",
// 24 : "Stopped",
// 28 : "VertexIsNotEndpointOfParent",

/* -- helpers -- */

dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec) {
  return dd4hep::rec::Vector3D(vec.x, vec.y, vec.z);
}
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3f vec) {
  return dd4hep::rec::Vector3D(static_cast<double>(vec.x), static_cast<double>(vec.y), static_cast<double>(vec.z));
}
edm4hep::Vector3d ConvertVector(dd4hep::rec::Vector3D vec) {
  return edm4hep::Vector3d(vec.x(), vec.y(), vec.z());
}

dd4hep::Position ConvertVector_toPosition(dd4hep::rec::Vector3D vec) {
  return dd4hep::Position(vec.x(), vec.y(), vec.z());
}
dd4hep::Position ConvertVector_toPosition(edm4hep::Vector3f vec) {
  return dd4hep::Position(vec.x, vec.y, vec.z);
}
dd4hep::Position ConvertVector_toPosition(edm4hep::Vector3d vec) {
  return dd4hep::Position(vec.x, vec.y, vec.z);
}

TGeoHMatrix ComputeSensorTrafoMatrix(const dd4hep::DDSegmentation::VolumeID& volumeID, const dd4hep::VolumeManager& volumeManager, const TGeoRotation& sensorNormalRotation) {
  TGeoHMatrix M = volumeManager.lookupDetElement(volumeID).nominal().worldTransformation();

  /* rotate the local coordinate system st. sensor U is (1,0,0), V is (0,1,0) and normal vector is (0,0,1) */
  M.Multiply(sensorNormalRotation);

  /* rotation is unitless, but need to convert translation from cm to mm (dd4hep::mm = 0.1) */
  double* transl = M.GetTranslation();
  transl[0] = transl[0] / dd4hep::mm;
  transl[1] = transl[1] / dd4hep::mm;
  transl[2] = transl[2] / dd4hep::mm;
  M.SetTranslation(transl);

  return M;
}

dd4hep::rec::Vector3D GlobalToLocal(const dd4hep::rec::Vector3D& global, const TGeoHMatrix& M) {
  double local[3];
  M.MasterToLocal(global, local);
  return dd4hep::rec::Vector3D(local[0], local[1], local[2]);
}

dd4hep::rec::Vector3D LocalToGlobal(const dd4hep::rec::Vector3D& local, const TGeoHMatrix& M) {
  double global[3];
  M.LocalToMaster(local, global);
  return dd4hep::rec::Vector3D(global[0], global[1], global[2]);
}

int GetLayer(const dd4hep::DDSegmentation::VolumeID& volumeID, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) {
  return static_cast<int>(cellIdDecoder->get(volumeID, "layer"));
}



/* -- Binning things -- */

int ComputeBinIndex(float x, float binX0, float binWidth, int binN) {
  /** Get the bin index for a given x value
   *  binX0 is the lower edge of the first bin
   *  binWidth is the width of the bins
   *  binN is the number of bins
   *  return -1 if x is out of range
   */

  if (binN <= 0) throw std::runtime_error("VTXdigi_tools::ComputeBinIndex(): binN must be positive");
  if (binWidth <= 0.0) throw std::runtime_error("VTXdigi_tools::ComputeBinIndex(): binWidth must be positive");

  float relativePos = (x - binX0) / binWidth; // shift to [0, binN]
  if (relativePos < 0.0f || relativePos > static_cast<float>(binN))
    return -1;
  if (relativePos == static_cast<float>(binN))
    return binN - 1; // include upper edge in last bin (makes sense for pixels)
  return static_cast<int>(relativePos);
} // ComputeBinIndex()

float ComputeBinCenter(int i, float binX0, float binWidth) {
  if (i < 0) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): bin index must be non-negative");
  if (binWidth <= 0.0) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): binWidth must be positive");

  return binX0 + (static_cast<float>(i) + 0.5f) * binWidth; // add 0.5*binWidth to shift from lower edge to center
} // ComputeBinCenter()
float ComputeBinCenter(int i, float binX0, float binX1, int binN) {
  if (binN <= 0) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): binN must be positive");
  if (binX1 <= binX0) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): binX1 must be greater than binX0");
  if (i < 0 || i >= binN) throw std::runtime_error("VTXdigi_tools::ComputeBinCenter(): bin index out of bounds");

  const float binWidth = (binX1 - binX0) / static_cast<float>(binN);
  return ComputeBinCenter(i, binX0, binWidth);
} // ComputeBinCenter()

std::pair<int, int> ComputePixelIndices(const dd4hep::rec::Vector3D& pos, const std::pair<float, float> pixelPitch, const std::pair<size_t, size_t> pixelCount) {
  const float length_u_half = 0.5 * pixelPitch.first * pixelCount.first;
  int i_u = ComputeBinIndex(
    pos.x(),
    -length_u_half,
    pixelPitch.first,
    pixelCount.first);

  const float length_v_half = 0.5 * pixelPitch.second * pixelCount.second;
  int i_v = ComputeBinIndex(
    pos.y(),
    -length_v_half,
    pixelPitch.second,
    pixelCount.second);

  return {i_u, i_v};
} // ComputePixelIndices()

std::array<int, 3> ComputeInPixelIndices(const dd4hep::rec::Vector3D& pos, const std::array<int, 3>& binCount, const std::pair<float, float>& pixelPitch, const std::array<float, 3>& activeVolumeDimensions) {
  std::array<int, 3> indices;

  const float posShifted_u = pos.x() + 0.5 * activeVolumeDimensions[0]; // shift to [0, length_u]
  if (posShifted_u < 0.0 || posShifted_u > activeVolumeDimensions[0]) {
    indices[0] = -1; // out of bounds
  }
  else {
    float posInPixel_u = std::fmod(posShifted_u,  pixelPitch.first);
    if (posInPixel_u < 0.0) posInPixel_u +=  pixelPitch.first; // ensure positive remainder
    indices[0] = ComputeBinIndex(posInPixel_u, 0.0,  pixelPitch.first / binCount[0], binCount[0]);
  }

  const float posShifted_v = pos.y() + 0.5 * activeVolumeDimensions[1];
  if (posShifted_v < 0.0 || posShifted_v > activeVolumeDimensions[1]) {
    indices[1] = -1; // out of bounds
  }
  else {
    float posInPixel_v = std::fmod(posShifted_v, pixelPitch.second);
    if (posInPixel_v < 0.0) posInPixel_v += pixelPitch.second;
    indices[1] = ComputeBinIndex(posInPixel_v, 0.0, pixelPitch.second / binCount[1], binCount[1]);
  }

  // vertical (w) binning: shift to [0, thickness]
  const float posShifted_w = pos.z() + 0.5 * activeVolumeDimensions[2];
  indices[2] = ComputeBinIndex(posShifted_w, 0.0, activeVolumeDimensions[2] / binCount[2], binCount[2]); // no fmod, so out-of-bounds is caught

  return indices;
} // ComputeInPixelIndices()

dd4hep::rec::Vector3D ComputePosFromPixIndex_local(const std::pair<int, int> pixelIndex, const std::pair<float, float> sensorLength,  const std::pair<float, float> pixelPitch, float depletedRegionDepthCenter) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */

  float u = (static_cast<float>(pixelIndex.first) + 0.5f) * pixelPitch.first - 0.5f * sensorLength.first; // in mm
  float v = (static_cast<float>(pixelIndex.second) + 0.5f) * pixelPitch.second - 0.5f * sensorLength.second;
  float w = depletedRegionDepthCenter;

  return dd4hep::rec::Vector3D(u, v, w);
}

dd4hep::rec::Vector3D ComputePosFromPixIndex_local(const std::pair<int, int> pixelIndex, const std::pair<float, float> sensorLength, const std::pair<float, float> pixelPitch) {
  return ComputePosFromPixIndex_local(pixelIndex, sensorLength, pixelPitch, 0.f);
}

dd4hep::rec::Vector3D ComputePosFromPixIndex_local(const std::pair<float, float> index, const std::pair<float, float> sensorLength,  const std::pair<float, float> pixelPitch, float depletedRegionDepthCenter) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */

  float u = (index.first + 0.5f) * pixelPitch.first - 0.5f * sensorLength.first; // in mm. Add 0.5*pixelPitch to shift from pixel edge to center, since index 0 is defined as the center of the pixel.
  float v = (index.second + 0.5f) * pixelPitch.second - 0.5f * sensorLength.second;
  float w = depletedRegionDepthCenter;

  return dd4hep::rec::Vector3D(u, v, w);
}

dd4hep::rec::Vector3D ComputePosFromPixIndex_local(const std::pair<float, float> index, const std::pair<float, float> sensorLength,  const std::pair<float, float> pixelPitch) {
  return ComputePosFromPixIndex_local(index, sensorLength, pixelPitch, 0.f);
}

/* -- HitMap -- */

HitMap::HitMap(std::pair<size_t, size_t> pixelCount) : m_pixCount(pixelCount) {
  const int inverseOccupancy = 2000; // assume occupancy, 5e-4 is quite conservative for Z-run
  m_pixels.reserve(pixelCount.first * pixelCount.second / inverseOccupancy); // avoid too many reallocations
}

void HitMap::FillCharge(std::pair<int, int> i_uv, float charge, const SimHitWrapper& simHitWrapper) {
  if (charge < 1.e-6f)
    return; // skip very small charge additions for performance (this is NECESSARY to skip in-pix bins with weight ~0)
  if (_OutOfBounds(i_uv)) [[unlikely]]
    throw std::runtime_error("HitMap::FillCharge: pixel i_u or i_v ( " + std::to_string(i_uv.first) + ", " + std::to_string(i_uv.second) + ") out of range");

  auto [iter, inserted] = m_pixels.try_emplace(i_uv, Pixel(i_uv));
  iter->second.charge += charge;
  iter->second.simHits.insert(&simHitWrapper);
}

void HitMap::ApplyChargeSmearing(const Rndm::Numbers& rndm_charge) {
  auto hitIter = m_pixels.begin();
  while (hitIter != m_pixels.end()) {
    hitIter->second.charge = std::max(hitIter->second.charge + static_cast<float>(rndm_charge()), 0.f); // don't allow negative charge after smearing
    ++hitIter;
  }
}

void HitMap::ApplyThreshold(const float threshold, const Rndm::Numbers* rndm_threshold) {
  auto hitIter = m_pixels.begin();
  while (hitIter != m_pixels.end()) {
    // optionally disperse the threshold per pixel (drawn per event per sensor per pixel)
    const float pixThreshold = rndm_threshold ? threshold + static_cast<float>((*rndm_threshold)()) : threshold;
    if (hitIter->second.charge < pixThreshold)
      hitIter = m_pixels.erase(hitIter); // erase returns the iterator to the next element, so this is safe to do while iterating
    else
      ++hitIter;
  }
}

float HitMap::GetCharge(std::pair<int, int> i_uv) const {
  if (_OutOfBounds(i_uv)) [[unlikely]] {
    throw std::runtime_error("HitMap::GetCharge: pixel i_u or i_v ( " + std::to_string(i_uv.first) + ", " + std::to_string(i_uv.second) + ") out of range");
  }
  auto it = m_pixels.find(i_uv);
  if (it == m_pixels.end())
    return 0.f; // if pixel not found, charge is 0
  return it->second.charge;
}

float HitMap::GetTotalCharge() const {
  float totalCharge = 0.f;
  for (const auto& [i_uv, pixHit] : m_pixels) {
    totalCharge += pixHit.charge;
  }
  return totalCharge;
}

inline bool HitMap::_OutOfBounds(std::pair<int, int> i_uv) const {
  return (
    i_uv.first < 0
    || i_uv.first >= static_cast<int>(m_pixCount.first)
    || i_uv.second < 0
    || i_uv.second >= static_cast<int>(m_pixCount.second)
  );
}

/* -- Clusterization -- */

std::pair<float, float> Cluster::ComputePos() const {
  std::pair<float, float> pos{0.f, 0.f};
  for (const Pixel* pix : pixels) {
    pos.first += pix->index.first * pix->charge;
    pos.second += pix->index.second * pix->charge;
  }
  pos.first /= charge;
  pos.second /= charge;
  return pos;
}

int Cluster::GetSize(const int axis) const {
  int min = std::numeric_limits<int>::max();
  int max = std::numeric_limits<int>::min();

  if (axis == 0) { // u
    for (const Pixel* pix : pixels) {
      min = std::min(min, pix->index.first);
      max = std::max(max, pix->index.first);
    }
  }
  else if (axis == 1) { // v
    for (const Pixel* pix : pixels) {
      min = std::min(min, pix->index.second);
      max = std::max(max, pix->index.second);
    }
  }
  else {
    throw std::runtime_error("Cluster::GetClusterSize: axis must be 0 (u) or 1 (v), got " + std::to_string(axis));
  }

  return max - min + 1; // +1 because of counting: if min=max, cluster size is 1, not 0
}

std::pair<float, float> Cluster::ComputePosUncertainty_ChargeWeighted(const std::pair<float, float>& clusterPos) const {
  float sig2_u=0.f, sig2_v=0.f;
  for (const Pixel* pix : pixels) {
    float du = (pix->index.first - clusterPos.first);
    float dv = (pix->index.second - clusterPos.second);
    sig2_u += pix->charge * du * du;
    sig2_v += pix->charge * dv * dv;
  }
  sig2_u /= charge;
  sig2_v /= charge;
  return {std::sqrt(sig2_u), std::sqrt(sig2_v)};
}

std::pair<float, float> Cluster::ComputePosUncertainty_ChargeWeighted() const {
  return ComputePosUncertainty_ChargeWeighted(ComputePos());
}

float Cluster::GetSeedPixelCharge() const {
  float maxCharge = 0;
  for (const auto pixel : pixels) {
    if (pixel->charge > maxCharge) {
      maxCharge = pixel->charge;
    }
  }
  return maxCharge;
}


std::array<std::pair<int, int>, 4> GetDirectNeighbors(const std::pair<int, int>& i_uv) {
  return {{
    {i_uv.first - 1, i_uv.second}, // left
    {i_uv.first + 1, i_uv.second}, // right
    {i_uv.first, i_uv.second - 1}, // down
    {i_uv.first, i_uv.second + 1}  // up
  }};
}
std::array<std::pair<int, int>, 8> GetNeighbors(const std::pair<int, int>& i_uv) {
  return {{
    {i_uv.first - 1, i_uv.second}, // left
    {i_uv.first - 1, i_uv.second + 1}, // upper left
    {i_uv.first, i_uv.second + 1}, // up
    {i_uv.first + 1, i_uv.second + 1}, // upper right
    {i_uv.first + 1, i_uv.second}, // right
    {i_uv.first + 1, i_uv.second - 1}, // lower right
    {i_uv.first, i_uv.second - 1}, // lower
    {i_uv.first - 1, i_uv.second - 1}, // lower left
  }};
}


std::vector<Cluster> HitMap::ComputeClusters_singePixels() const {
  std::vector<Cluster> clusters;

  for (const auto& p : m_pixels) {
    const Pixel* pixel = &(p.second); // cluster stores pointers to pixels (pixels are stored in HitMap as values)

    clusters.emplace_back();
    clusters.back().pixels.push_back(pixel);
    clusters.back().charge = pixel->charge;
    for (const SimHitWrapper* simHitWrapper : pixel->simHits) {
      clusters.back().simHits.insert(simHitWrapper);
    }
  } // loop over pixelHits

  return clusters;
}

std::vector<Cluster> HitMap::ComputeClusters() const {
  /* Breadth First Search (BFS) implementation for clustering */

  std::vector<Cluster> clusters;
  std::unordered_set<std::pair<int,int>, Hash_PairInt> visited;

  for (const auto& p : m_pixels) {
    const std::pair<int,int> seed_uv = p.first;
    if (visited.contains(seed_uv))
      continue;

    clusters.emplace_back(); // create new cluster
    clusters.back().pixels.reserve(10); // 10 should include >90% of clusters. i guess.

    std::queue<std::pair<int,int>> queue;
    queue.push(seed_uv);
    visited.insert(seed_uv);

    while (!queue.empty()) {
      const std::pair<int,int> current_uv = queue.front();
      queue.pop();

      /* Add pixl to cluster */
      const Pixel* pixel = &(m_pixels.at(current_uv)); // get pixel pointer from map
      clusters.back().pixels.push_back(pixel);
      clusters.back().charge += pixel->charge;
      for (const SimHitWrapper* simHitWrapper : pixel->simHits) {
        clusters.back().simHits.insert(simHitWrapper);
      }
      /* Add all neighboring pixels to queue */
      // for (const auto& neighbor_uv : GetDirectNeighbors(current_uv)) {
      for (const auto& neighbor_uv : GetNeighbors(current_uv)) {
        if (!m_pixels.contains(neighbor_uv))
          continue;
        if (visited.contains(neighbor_uv))
          continue;
        queue.push(neighbor_uv);
        visited.insert(neighbor_uv);
      }

    } // loop over queue
  } // loop over cluster-seeds
  return clusters;
}

} // namespace VTXdigi_tools
