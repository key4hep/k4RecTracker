// VTXdigi_Modular/src/VTXdigi_tools.cpp
#include "VTXdigi_tools.h"
#include "VTXdigi_Modular.h"

#include <DD4hep/Objects.h>
#include <DD4hep/VolumeManager.h>
#include <DDRec/Vector3D.h>
#include <DDSegmentation/Segmentation.h>
#include <Parsers/Primitives.h>
#include <array>
#include <edm4hep/Vector3d.h>

namespace VTXdigi_tools {

SimHitWrapper::SimHitWrapper(
  edm4hep::SimTrackerHit simTrackerHit, dd4hep::DDSegmentation::VolumeID volumeID, const VTXdigi_Modular& digitizer)
    : m_simTrackerHit(simTrackerHit), m_volumeID(volumeID) {

  m_charge = static_cast<float>(m_simTrackerHit.getEDep() * (dd4hep::GeV / dd4hep::keV) * kChargePerkeV); // convert energy deposit (in keV) to number of electrons
  m_layerNumber = digitizer.GetLayer(m_volumeID);
  // m_truthPos is set later in VTXdigi_Modular::operator() via SetTruthPos(...) to avoid double-calculating the sensor transformation matrix

  m_mcParticleLevel = ComputeMCParticleLevel(m_simTrackerHit, m_volumeID, digitizer);
}

void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept {
  std::swap(a.m_simTrackerHit, b.m_simTrackerHit);
  std::swap(a.m_volumeID, b.m_volumeID);
  std::swap(a.m_charge, b.m_charge);
  std::swap(a.m_layerNumber, b.m_layerNumber);
  std::swap(a.m_truthPos, b.m_truthPos);
  std::swap(a.m_mcParticleLevel, b.m_mcParticleLevel);
} // swap(Hit&, Hit&)

MCParticleLevel ComputeMCParticleLevel(const edm4hep::SimTrackerHit& simTrackerHit, dd4hep::DDSegmentation::VolumeID volumeID, const VTXdigi_Modular& digitizer) {
  if ( simTrackerHit.isProducedBySecondary() ) {
    // ddsim drops MCParticles below a certain energy cut to save computing cost and disk space.
    // so if ddsim dropped the MCParticle that caused this simHit, we assume it was a delta ray
    return MCParticleLevel::Delta;
  }
  else {
    const int32_t simulatorStatus = simTrackerHit.getParticle().getSimulatorStatus();
    const int32_t mask = 1 << edm4hep::MCParticle::BITCreatedInSimulation; // should be bit 30
    const bool causedByPrimary = (simulatorStatus & mask) == 0; // bit is not set -> created in generator
    if ( causedByPrimary ) {
      return MCParticleLevel::Primary;
    }
    else {
      // now check if the MCParticle prod. vertex lies outside this sensors volume (by comparing volumeIDs)
      const dd4hep::rec::Vector3D prodVertex = ConvertVector(simTrackerHit.getParticle().getVertex());
      const dd4hep::DDSegmentation::CellID prodVertex_cellID = digitizer.GetCellID(prodVertex);
      const dd4hep::DDSegmentation::VolumeID prodVertex_volumeID = digitizer.GetVolumeID(prodVertex_cellID);
      // VolumeID is 0 if pos is outside of any sensitive volume

      if (prodVertex_volumeID == 0 || prodVertex_volumeID != volumeID)
        return MCParticleLevel::Secondary; // MCParticle created outside of this sensor's sensitive volume
      else
        return MCParticleLevel::Delta;
    }
  }
}

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

std::string VectorToString(const dd4hep::rec::Vector3D& vec) {
  std::ostringstream oss;
  oss << "(" << vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
  return oss.str();
}

bool IsInsideVolume(const dd4hep::rec::Vector3D& pos, const std::array<double, 3>& dims, const double tolerance) {
  return std::abs(pos.x()) <= 0.5 * dims[0] + tolerance
      && std::abs(pos.y()) <= 0.5 * dims[1] + tolerance
      && std::abs(pos.z()) <= 0.5 * dims[2] + tolerance;
}

dd4hep::rec::Vector3D ClampToVolume(const dd4hep::rec::Vector3D& pos, const std::array<double, 3>& dims) {
  return dd4hep::rec::Vector3D(
    std::clamp(pos.x(), -0.5 * dims[0], 0.5 * dims[0]),
    std::clamp(pos.y(), -0.5 * dims[1], 0.5 * dims[1]),
    std::clamp(pos.z(), -0.5 * dims[2], 0.5 * dims[2]));
}

int GetLayer(const dd4hep::DDSegmentation::VolumeID& volumeID, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) {
  return static_cast<int>(cellIdDecoder->get(volumeID, "layer"));
}

dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec) {
  return dd4hep::rec::Vector3D(vec.x, vec.y, vec.z);
}
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3f vec) {
  return dd4hep::rec::Vector3D(static_cast<double>(vec.x), static_cast<double>(vec.y), static_cast<double>(vec.z));
}
edm4hep::Vector3d ConvertVector(dd4hep::rec::Vector3D vec) {
  return edm4hep::Vector3d(vec.x(), vec.y(), vec.z());
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

dd4hep::rec::Vector3D TrafoVec_global_local(const dd4hep::rec::Vector3D& global, const TGeoHMatrix& M) {
  double local[3];
  M.MasterToLocalVect(global, local);
  return dd4hep::rec::Vector3D(local[0], local[1], local[2]);
}

dd4hep::rec::Vector3D TrafoVec_local_global(const dd4hep::rec::Vector3D& local, const TGeoHMatrix& M) {
  double global[3];
  M.LocalToMasterVect(local, global);
  return dd4hep::rec::Vector3D(global[0], global[1], global[2]);
}



dd4hep::rec::Vector3D Trafo_global_local(const dd4hep::rec::Vector3D& global, const TGeoHMatrix& M) {
  double local[3];
  M.MasterToLocal(global, local);
  return dd4hep::rec::Vector3D(local[0], local[1], local[2]);
}

dd4hep::rec::Vector3D Trafo_local_global(const dd4hep::rec::Vector3D& local, const TGeoHMatrix& M) {
  double global[3];
  M.LocalToMaster(local, global);
  return dd4hep::rec::Vector3D(global[0], global[1], global[2]);
}

PixelCoords Trafo_local_pixCoords(const dd4hep::rec::Vector3D& local, const std::array<double, 2> pixelPitch, const std::array<size_t, 2> pixelCount) {
  const std::array<double, 2> local_2d = {local.x(), local.y()};
  PixelCoords pixC;
  for (size_t axis = 0; axis < 2; ++axis) {
    const double halfLength = 0.5 * pixelPitch[axis] * pixelCount[axis];
    const double clamped = std::clamp(local_2d[axis], -halfLength, halfLength);
    pixC[axis] = (clamped + halfLength) / pixelPitch[axis] - 0.5; // shift from [-halfLength, halfLength] to [-0.5, pixelCount - 0.5]
  }
  return pixC;
}

dd4hep::rec::Vector3D Trafo_pixCoords_local(const PixelCoords pixCoords,  const std::array<double, 2> pixelPitch, const std::array<size_t, 2> pixelCount, const double w) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */
  double u = (pixCoords[0] + 0.5) * pixelPitch[0] - 0.5 * pixelPitch[0] * pixelCount[0]; // in mm. Add 0.5*pixelPitch to shift from pixel edge to center, since index 0 is defined as the center of the pixel.
  double v = (pixCoords[1] + 0.5) * pixelPitch[1] - 0.5 * pixelPitch[1] * pixelCount[1];

  return dd4hep::rec::Vector3D(u, v, w);
}


/* -- Binning things -- */

int ComputeBinIndex(double x, double binX0, double binWidth, int binN) {
  #ifndef NDEBUG
    if (binN <= 0) throw std::runtime_error("VTXdigi_tools::ComputeBinIndex(): binN must be positive");
    if (binWidth <= 0.0) throw std::runtime_error("VTXdigi_tools::ComputeBinIndex(): binWidth must be positive");
  #endif

  const double relativePos = (x - binX0) / binWidth; // shift to [0, binN]
  return std::clamp(static_cast<int>(std::floor(relativePos)), 0, binN - 1);
} // ComputeBinIndex()

PixelIndex Trafo_local_pixI(const dd4hep::rec::Vector3D& local, const std::array<double, 2> pixelPitch, const std::array<size_t, 2> pixelCount) {
  PixelIndex pixI;
  const double length_u_half = 0.5 * pixelPitch[0] * pixelCount[0];
  pixI[0] = ComputeBinIndex(
    local.x(),
    -length_u_half,
    pixelPitch[0],
    pixelCount[0]);

  const double length_v_half = 0.5 * pixelPitch[1] * pixelCount[1];
  pixI[1] = ComputeBinIndex(
    local.y(),
    -length_v_half,
    pixelPitch[1],
    pixelCount[1]);

  return pixI;
} // Trafo_local_pixI()

std::pair<PixelIndex, VoxelIndex> Trafo_local_pixIVoxI(const dd4hep::rec::Vector3D& local, const std::array<double, 2> pixelPitch, const std::array<size_t, 2> pixelCount, const double sensorActiveThickness, const std::array<int, 3> voxelCount) {
  const std::array<double, 2> local_uv = {local.x(), local.y()};
  PixelIndex pixI;
  VoxelIndex voxI;

  // binning in u/v: pixel + voxel grid
  for (int axis = 0; axis < 2; ++axis) {
    const int g = ComputeBinIndex(
      local_uv[axis],
      -0.5 * pixelPitch[axis] * pixelCount[axis],
      pixelPitch[axis] / voxelCount[axis],
      static_cast<int>(pixelCount[axis]) * voxelCount[axis]);
    pixI[axis] = g / voxelCount[axis];
    voxI[axis] = g % voxelCount[axis];
  }

  // binning in w: only voxel grid
  voxI[2] = ComputeBinIndex(local.z(), -0.5 * sensorActiveThickness, sensorActiveThickness / voxelCount[2], voxelCount[2]);

  return {pixI, voxI};
} // Trafo_local_pixIVoxI()

dd4hep::rec::Vector3D Trafo_pixI_local(const PixelIndex pixI, const std::array<double, 2> pixelPitch, const std::array<size_t, 2> pixelCount, const double w) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */

  double u = (pixI[0] + 0.5) * pixelPitch[0] - 0.5 * pixelPitch[0] * pixelCount[0]; // in mm
  double v = (pixI[1] + 0.5) * pixelPitch[1] - 0.5 * pixelPitch[1] * pixelCount[1];

  return dd4hep::rec::Vector3D(u, v, w);
}

/* -- HitMap -- */

HitMap::HitMap(std::array<size_t, 2> pixelCount) : m_pixCount(pixelCount) {
  const int inverseOccupancy = 2000; // assume occupancy, 5e-4 is quite conservative for Z-run
  m_pixels.reserve(pixelCount[0] * pixelCount[1] / inverseOccupancy); // avoid too many reallocations
}

void HitMap::FillCharge(PixelIndex pixI, float charge, const SimHitWrapper& simHitWrapper) {
  if (charge < 1.e-6f)
    return; // skip very small charge additions for performance (this is NECESSARY to skip in-pix bins with weight ~0)
  if (_OutOfBounds(pixI)) [[unlikely]]
    throw std::runtime_error("HitMap::FillCharge: pixel i_u or i_v ( " + std::to_string(pixI[0]) + ", " + std::to_string(pixI[1]) + ") out of range");

  auto [iter, inserted] = m_pixels.try_emplace(pixI, Pixel(pixI));
  iter->second.charge += charge;
  iter->second.simHits.insert(&simHitWrapper);
}

void HitMap::ApplyChargeSmearing(const float sigma, TRandom3& randomGen) {
  auto hitIter = m_pixels.begin();
  while (hitIter != m_pixels.end()) {
    hitIter->second.charge = std::max(hitIter->second.charge + static_cast<float>(randomGen.Gaus(0, sigma)), 0.f); // don't allow negative charge after smearing
    ++hitIter;
  }
}

void HitMap::ApplyThreshold(const float threshold, const float thresholdDispersion, TRandom3& randomGen) {
  auto pixelIter = m_pixels.begin();
  while (pixelIter != m_pixels.end()) {
    // skip expensive Mersenne Twister step for pixels that are well above threshold
    if (pixelIter->second.charge > threshold + 8.f * thresholdDispersion) {
      ++pixelIter;
      continue;
    }

    float pixThreshold = threshold;
    if (thresholdDispersion != 0.f)
      pixThreshold += static_cast<float>(randomGen.Gaus(0, thresholdDispersion));

    if (pixelIter->second.charge < pixThreshold)
      pixelIter = m_pixels.erase(pixelIter); // erase returns the iterator to the next element
    else
      ++pixelIter;
  }
}


float HitMap::GetCharge(PixelIndex pixI) const {
  if (_OutOfBounds(pixI)) [[unlikely]] {
    throw std::runtime_error("HitMap::GetCharge: pixel i_u or i_v ( " + std::to_string(pixI[0]) + ", " + std::to_string(pixI[1]) + ") out of range");
  }
  auto it = m_pixels.find(pixI);
  if (it == m_pixels.end())
    return 0.f; // if pixel not found, charge is 0
  return it->second.charge;
}

float HitMap::GetTotalCharge() const {
  float totalCharge = 0.f;
  for (const auto& [pixI, pixHit] : m_pixels) {
    totalCharge += pixHit.charge;
  }
  return totalCharge;
}

inline bool HitMap::_OutOfBounds(PixelIndex pixI) const {
  return (
    pixI[0] < 0
    || pixI[0] >= static_cast<int>(m_pixCount[0])
    || pixI[1] < 0
    || pixI[1] >= static_cast<int>(m_pixCount[1])
  );
}

/* -- Clusterization -- */

PixelCoords Cluster::ComputeCoG(const bool clusterizeEndPixelsOnly) const {
  if (pixels.empty())
    throw std::runtime_error("Cluster::ComputeCoG: cluster has no pixels");

  PixelCoords pixCoords{0.f, 0.f};
  if (!clusterizeEndPixelsOnly) {
    for (const Pixel* pix : pixels) {
      pixCoords[0] += pix->index[0] * pix->charge;
      pixCoords[1] += pix->index[1] * pix->charge;
    }
    pixCoords[0] /= charge;
    pixCoords[1] /= charge;
  }
  else {
    for (int axis = 0; axis < 2; ++axis) {
      // find pixels with min and max index along the axis
      int index_min = std::numeric_limits<int>::max();
      int index_max = std::numeric_limits<int>::min();
      for (const Pixel* pix : pixels) {
        index_min = std::min(index_min, pix->index[axis]);
        index_max = std::max(index_max, pix->index[axis]);
      }
      // compute charge-weighted average of the min and max pixels
      float charge_min=0.f, charge_max=0.f;
      for (const Pixel* pix : pixels) {
        if (pix->index[axis] == index_min) charge_min += pix->charge;
        if (pix->index[axis] == index_max) charge_max += pix->charge;
      }
      float offset = (charge_max - charge_min) / (charge_min + charge_max) / 2.f; // offset in range [-0.5, 0.5] to shift the CoG towards the pixel with more charge

      pixCoords[axis] = static_cast<float>(index_min + index_max) * 0.5f + offset;
    }
  }
  return pixCoords;
}

int Cluster::GetSize(const int axis) const {
  int min = std::numeric_limits<int>::max();
  int max = std::numeric_limits<int>::min();

  if (axis == 0) { // u
    for (const Pixel* pix : pixels) {
      min = std::min(min, pix->index[0]);
      max = std::max(max, pix->index[0]);
    }
  }
  else if (axis == 1) { // v
    for (const Pixel* pix : pixels) {
      min = std::min(min, pix->index[1]);
      max = std::max(max, pix->index[1]);
    }
  }
  else {
    throw std::runtime_error("Cluster::GetClusterSize: axis must be 0 (u) or 1 (v), got " + std::to_string(axis));
  }

  return max - min + 1; // +1 because of counting: if min=max, cluster size is 1, not 0
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


std::array<PixelIndex, 4> GetDirectNeighbors(const PixelIndex& pixI) {
  return {{
    {pixI[0] - 1, pixI[1]}, // left
    {pixI[0] + 1, pixI[1]}, // right
    {pixI[0], pixI[1] - 1}, // down
    {pixI[0], pixI[1] + 1}  // up
  }};
}
std::array<PixelIndex, 8> GetNeighbors(const PixelIndex& pixI) {
  return {{
    {pixI[0] - 1, pixI[1]}, // left
    {pixI[0] - 1, pixI[1] + 1}, // upper left
    {pixI[0], pixI[1] + 1}, // up
    {pixI[0] + 1, pixI[1] + 1}, // upper right
    {pixI[0] + 1, pixI[1]}, // right
    {pixI[0] + 1, pixI[1] - 1}, // lower right
    {pixI[0], pixI[1] - 1}, // lower
    {pixI[0] - 1, pixI[1] - 1}, // lower left
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
  std::unordered_set<PixelIndex, Hash_PixIndex> visited;

  for (const auto& p : m_pixels) {
    const PixelIndex seedI = p.first;
    if (visited.contains(seedI))
      continue;

    clusters.emplace_back(); // create new cluster
    clusters.back().pixels.reserve(10); // 10 should include >90% of clusters. i guess.

    std::queue<PixelIndex> queue;
    queue.push(seedI);
    visited.insert(seedI);

    while (!queue.empty()) {
      const PixelIndex currentI = queue.front();
      queue.pop();

      /* Add pixl to cluster */
      const Pixel* pixel = &(m_pixels.at(currentI)); // get pixel pointer from map
      clusters.back().pixels.push_back(pixel);
      clusters.back().charge += pixel->charge;
      for (const SimHitWrapper* simHitWrapper : pixel->simHits) {
        clusters.back().simHits.insert(simHitWrapper);
      }
      /* Add all neighboring pixels to queue */
      for (const auto& neighborI : GetNeighbors(currentI)) {
        if (!m_pixels.contains(neighborI))
          continue;
        if (visited.contains(neighborI))
          continue;
        queue.push(neighborI);
        visited.insert(neighborI);
      }
    } // loop over queue
  } // loop over cluster-seeds
  return clusters;
}

} // namespace VTXdigi_tools
