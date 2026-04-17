// VTXdigi_Modular/src/VTXdigi_tools.cpp
#include "VTXdigi_tools.h"

namespace VTXdigi_tools {

SimHitWrapper::SimHitWrapper(edm4hep::SimTrackerHit simTrackerHit, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) : m_simTrackerHit(simTrackerHit) {
  m_cellID = GetCellID_short(m_simTrackerHit);

  m_charge = static_cast<float>(m_simTrackerHit.getEDep() * (dd4hep::GeV / dd4hep::keV) * kChargePerkeV); // convert energy deposit (in keV) to number of electrons 

  m_layerNumber = GetLayer(m_cellID, cellIdDecoder);
}

void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept {
  std::swap(a.m_simTrackerHit, b.m_simTrackerHit);
  std::swap(a.m_cellID, b.m_cellID);
  std::swap(a.m_charge, b.m_charge);
  std::swap(a.m_layerNumber, b.m_layerNumber);
} // swap(Hit&, Hit&)

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

TGeoHMatrix ComputeSensorTrafoMatrix(const dd4hep::DDSegmentation::CellID& cellID, const dd4hep::VolumeManager& volumeManager, const TGeoRotation& sensorNormalRotation) {
  TGeoHMatrix M = volumeManager.lookupDetElement(cellID).nominal().worldTransformation();

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

dd4hep::DDSegmentation::CellID GetCellID_short(const edm4hep::SimTrackerHit& simTrackerHit) {
  /* Mask removes segmentation bits, now cellID is unique for each sensor */
  std::uint64_t m_mask = (static_cast<std::uint64_t>(1) << 32) - 1;
  return simTrackerHit.getCellID() & m_mask;
}
dd4hep::DDSegmentation::CellID GetCellID_short(const dd4hep::DDSegmentation::CellID& cellID) {
  /* Mask removes segmentation bits, now cellID is unique for each sensor */
  std::uint64_t m_mask = (static_cast<std::uint64_t>(1) << 32) - 1;
  return cellID & m_mask;
}

int GetLayer(const dd4hep::DDSegmentation::CellID& cellID, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) {
  return static_cast<int>(cellIdDecoder->get(cellID, "layer"));
} 
int GetLayer(const edm4hep::SimTrackerHit& simTrackerHit, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) {
  return GetLayer(GetCellID_short(simTrackerHit), cellIdDecoder);
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

std::array<int, 3> ComputeInPixelIndices(const dd4hep::rec::Vector3D& pos, const std::array<int, 3>& binCount, const std::pair<float, float>& pixelPitch, const std::array<float, 3>& sensorDimensions) {
  std::array<int, 3> indices;

  const float posShifted_u = pos.x() + 0.5 * sensorDimensions[0]; // shift to [0, length_u]
  if (posShifted_u < 0.0 || posShifted_u > sensorDimensions[0]) {
    indices[0] = -1; // out of bounds
  }
  else {
    float posInPixel_u = std::fmod(posShifted_u,  pixelPitch.first);
    if (posInPixel_u < 0.0) posInPixel_u +=  pixelPitch.first; // ensure positive remainder
    indices[0] = ComputeBinIndex(posInPixel_u, 0.0,  pixelPitch.first / binCount[0], binCount[0]);
  }

  const float posShifted_v = pos.y() + 0.5 * sensorDimensions[1];
  if (posShifted_v < 0.0 || posShifted_v > sensorDimensions[1]) {
    indices[1] = -1; // out of bounds
  }
  else {
    float posInPixel_v = std::fmod(posShifted_v, pixelPitch.second);
    if (posInPixel_v < 0.0) posInPixel_v += pixelPitch.second;
    indices[1] = ComputeBinIndex(posInPixel_v, 0.0, pixelPitch.second / binCount[1], binCount[1]);
  }

  // vertical (w) binning: shift to [0, thickness]
  const float posShifted_w = pos.z() + 0.5 * sensorDimensions[2];
  indices[2] = ComputeBinIndex(posShifted_w, 0.0, sensorDimensions[2] / binCount[2], binCount[2]); // no fmod, so out-of-bounds is caught

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

void HitMap::ApplyThreshold(const float threshold) {
  auto hitIter = m_pixels.begin();
  while (hitIter != m_pixels.end()) {
    if (hitIter->second.charge < threshold)
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


std::array<std::pair<int, int>, 4> GetDirectNeighbors(const std::pair<int, int>& i_uv) {
  return {{
    {i_uv.first - 1, i_uv.second}, // left
    {i_uv.first + 1, i_uv.second}, // right
    {i_uv.first, i_uv.second - 1}, // down
    {i_uv.first, i_uv.second + 1}  // up
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
      for (const auto& neighbor_uv : GetDirectNeighbors(current_uv)) {
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

/* -- Tool tests -- */

bool ToolTest() {
  std::cout << " | Running VTXdigi tool tests" << std::endl;
  bool passed = true;
  
  std::cout << " | VTXdigi_tools::ComputeBinIndex()";
  {
    bool passedInternal = true;
    const float binX0 = -1.0f;
    const float binWidth = 0.5f;
    const int binN = 6;

    const float inputs[5] = { -1.0f, 0.0f, 2.0f, -1.1f, 2.1f };
    const int expectedOutputs[5] = { 0, 2, 5, -1, -1 };

    for (size_t i = 0; i < 5; ++i) {
      int result = ComputeBinIndex(inputs[i], binX0, binWidth, binN);
      if (result != expectedOutputs[i]) {
        std::cout << " - FAILED " << std::endl << " | -> Expected bin index " << expectedOutputs[i] << " for x=" << inputs[i] << ", got " << result << std::endl;
        passedInternal = false;
      }
    }
    if (passedInternal)
      std::cout << " - PASSED" << std::endl;
    passed = passed && passedInternal;
  }

  std::cout << " | VTXdigi_tools::ComputePixelIndices()";
  {
    bool passedInternal = true;

    const std::pair<float, float> pixelPitch = { 1.0f, 2.0f };
    const std::pair<size_t, size_t> pixelCount = { 10, 10 };

    std::array<dd4hep::rec::Vector3D, 10> inputs = {{dd4hep::rec::Vector3D(0.f, 0.f, 0.f)}};
    inputs = {
      dd4hep::rec::Vector3D( -4.5f, -9.0f, 0.0f ),
      dd4hep::rec::Vector3D(0.0f, 0.0f, 0.0f),
      dd4hep::rec::Vector3D(4.4f, 8.9f, 0.0f),
      dd4hep::rec::Vector3D(-5.0f, -10.0f, 0.0f), // edge tests
      dd4hep::rec::Vector3D(5.0f, 10.0f, 0.0f),
      dd4hep::rec::Vector3D(-5.1f, 0.0f, 0.0f), // out of bounds tests
      dd4hep::rec::Vector3D(5.1f, 0.0f, 0.0f), 
      dd4hep::rec::Vector3D(0.0f, -10.1f, 0.0f),
      dd4hep::rec::Vector3D(0.0f, 10.1f, 0.0f),
      dd4hep::rec::Vector3D(5.1f, 10.1f, 0.0f),
    };
    std::array<std::pair<int, int>, inputs.size()> expectedOutputs = {{{0, 0}}};
    expectedOutputs= {{
      {0,0},
      {5,5},
      {9,9},
      {0,0},
      {9,9},
      {-1,5},
      {-1,5},
      {5,-1},
      {5,-1},
      {-1,-1}
    }};

    for (size_t i = 0; i < inputs.size(); ++i) {
      std::pair<int, int> result = ComputePixelIndices(inputs.at(i), pixelPitch, pixelCount);
      if (result != expectedOutputs[i]) {
        std::cout << " - FAILED " << std::endl << " | -> Expected pixel indices (" << expectedOutputs[i].first << ", " << expectedOutputs[i].second << ") for (u,v,w)=(" << inputs.at(i).x() << ", " << inputs.at(i).y() << ", " << inputs.at(i).z() << "), got (" << result.first << ", " << result.second << ")" << std::endl;
        passedInternal = false;
      }
    }

    if (passedInternal)
      std::cout << " - PASSED" << std::endl;
    passed = passed && passedInternal;
  }

  std::cout << " | VTXdigi_tools::ComputePosFromPixIndex_local(std::pair<int,int>)";
  {
    bool passedInternal = true;

    const std::pair<float, float> sensorLength = { 10.0, 20.0 }; // 10 x 10 pixels
    const std::pair<float, float> pixelPitch = { 1.0, 2.0 };

    std::array<std::pair<int, int>, 6> inputs = {{{0, 0}}};
    inputs = {{
      {0, 0},
      {4, 4},
      {5, 5},
      {9, 9},
      {-1, -1}, // out of bounds test
      {10, 10},
    }};
    std::array<dd4hep::rec::Vector3D, inputs.size()> expectedOutputs = {{dd4hep::rec::Vector3D(0.f, 0.f, 0.f)}};
    expectedOutputs = {{
      dd4hep::rec::Vector3D( -4.5, -9.0, 0.0 ),
      dd4hep::rec::Vector3D( -0.5, -1.0, 0.0 ),
      dd4hep::rec::Vector3D( 0.5, 1.0, 0.0 ),
      dd4hep::rec::Vector3D( 4.5, 9.0, 0.0 ),
      dd4hep::rec::Vector3D( -5.5, -11.0, 0.0 ),
      dd4hep::rec::Vector3D( 5.5, 11.0, 0.0 ),
    }};
    
    for (size_t i = 0; i < inputs.size(); ++i) {
      dd4hep::rec::Vector3D result = ComputePosFromPixIndex_local(inputs.at(i), sensorLength, pixelPitch);
      if (std::abs(result.x() - expectedOutputs[i].x()) > 1e-6 || std::abs(result.y() - expectedOutputs[i].y()) > 1e-6) {
        if (passedInternal)
        std::cout << " - FAILED " << std::endl;
        std::cout << " | -> Expected local position (" << expectedOutputs[i].x() << ", " << expectedOutputs[i].y() << ") for index (i_u,i_v)=(" << inputs.at(i).first << ", " << inputs.at(i).second << "), got (" << result.x() << ", " << result.y() << ")" << std::endl;
        passedInternal = false;
      }
    }

    if (passedInternal)
      std::cout << " - PASSED" << std::endl;
    passed = passed && passedInternal;
  }

  std::cout << " | VTXdigi_tools::ComputePosFromPixIndex_local(std::pair<float,float>)";
  {
    bool passedInternal = true;

    const std::pair<float, float> pixelPitch = { 1.0, 2.0 };
    const std::pair<float, float> sensorLength = { 10.0, 20.0 }; // 10 x 10 pixels

    std::array<std::pair<float, float>, 9> inputs = {{{0., 0.}}};
    inputs = {{
      {0., 0.},
      {4., 4.},
      {4.5, 4.5},
      {5., 5.},
      {9.0, 9.0},
      {-0.5, -0.5},
      {9.5, 9.5},
      {-1., -1.}, // expect out of bounds to not be treated at all (done elsewhere)
      {10., 10.},
    }};
    std::array<dd4hep::rec::Vector3D, inputs.size()> expectedOutputs = {{dd4hep::rec::Vector3D(0.f, 0.f, 0.f)}};
    expectedOutputs = {{
      dd4hep::rec::Vector3D( -4.5, -9.0, 0.0 ),
      dd4hep::rec::Vector3D( -0.5, -1.0, 0.0 ),
      dd4hep::rec::Vector3D( 0., 0., 0. ),
      dd4hep::rec::Vector3D( 0.5, 1.0, 0.0 ),
      dd4hep::rec::Vector3D( 4.5, 9.0, 0.0 ),
      dd4hep::rec::Vector3D( -5.0, -10.0, 0.0 ),
      dd4hep::rec::Vector3D( 5.0, 10.0, 0.0 ),
      dd4hep::rec::Vector3D( -5.5, -11.0, 0.0 ),
      dd4hep::rec::Vector3D( 5.5, 11.0, 0.0 ),
    }};
    
    for (size_t i = 0; i < inputs.size(); ++i) {
      dd4hep::rec::Vector3D result = ComputePosFromPixIndex_local(inputs.at(i), sensorLength, pixelPitch);
      if (std::abs(result.x() - expectedOutputs[i].x()) > 1e-6 || std::abs(result.y() - expectedOutputs[i].y()) > 1e-6) {
        if (passedInternal)
        std::cout << " - FAILED " << std::endl;
        std::cout << " | -> Expected local position (" << expectedOutputs[i].x() << ", " << expectedOutputs[i].y() << ") for index (i_u,i_v)=(" << inputs.at(i).first << ", " << inputs.at(i).second << "), got (" << result.x() << ", " << result.y() << ")" << std::endl;
        passedInternal = false;
      }
    }

    if (passedInternal)
      std::cout << " - PASSED" << std::endl;
    passed = passed && passedInternal;
  }

  std::cout << " | VTXdigi_tools::ComputeInPixelIndices()";
  {
    bool passedInternal = true;

    const std::array<int, 3> binCount = { 4, 4, 4 };
    const std::pair<float, float> pixelPitch = { 1.0f, 2.0f }; // bin-width is 0.25 

    const std::array<float, 3> sensorDimensions = { 10.0f, 20.0f, 1.0f };

    std::array<dd4hep::rec::Vector3D, 10> inputs = {{dd4hep::rec::Vector3D(0.f, 0.f, 0.f)}};
    inputs = {
      dd4hep::rec::Vector3D(0.1, 0.1, -0.4),
      dd4hep::rec::Vector3D(0.1, 1.9, 0.1),
      dd4hep::rec::Vector3D(-0.1, -0.1, -0.4),
      dd4hep::rec::Vector3D(-0.1, -0.1, 0.1),
      dd4hep::rec::Vector3D(0., 0., 0.), // edges
      dd4hep::rec::Vector3D(1., 2.,0.5),
      dd4hep::rec::Vector3D(-5.1, -10.1, 0.), // out-of-bounds
      dd4hep::rec::Vector3D(5.1, 10.1, 0.),
      dd4hep::rec::Vector3D(0.1, 0.1, -0.6),
      dd4hep::rec::Vector3D(0.1, 0.1, 0.6),
    };
    std::array<std::array<int, 3>, inputs.size()> expectedOutputs = {{{0, 0, 0}}};
    expectedOutputs = {{
      {0, 0, 0},
      {0, 3, 2},
      {3, 3, 0},
      {3, 3, 2},
      {0, 0, 2}, // edges
      {0, 0, 3},
      {-1, -1, 2}, // out-of-bounds
      {-1, -1, 2},
      {0, 0, -1},
      {0, 0, -1},
    }};

    for (size_t i = 0; i < inputs.size(); ++i) {
      std::array<int, 3> result = ComputeInPixelIndices(inputs.at(i), binCount, pixelPitch, sensorDimensions);
      if (result != expectedOutputs[i]) {
        std::cout << " - FAILED " << std::endl << " | -> Expected in-pixel indices (" << expectedOutputs[i][0] << ", " << expectedOutputs[i][1] << ", " << expectedOutputs[i][2] << ") for (u,v,w)=(" << inputs.at(i).x() << ", " << inputs.at(i).y() << ", " << inputs.at(i).z() << "), got (" << result[0] << ", " << result[1] << ", " << result[2] << ")" << std::endl;
        passedInternal = false;
      }
    }

    if (passedInternal)
      std::cout << " - PASSED" << std::endl;
    passed = passed && passedInternal;
  }

  std::cout << " | VTXdigi_tools::HitMap()";
  {
    bool passedInternal = true;

    SimHitWrapper simHitWrapper; // not used, just need something to pass to FillCharge()

    HitMap hitMap({10,10});

    hitMap.FillCharge({0,0},1, simHitWrapper);
    hitMap.FillCharge({1,1},2, simHitWrapper);

    if (hitMap.GetCharge({0,0}) != 1) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 1 at (0,0), got " << hitMap.GetCharge({0,0}) << std::endl;
      passedInternal = false;
    }
    if (hitMap.GetCharge({1,1}) != 2) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 2 at (1,1), got " << hitMap.GetCharge({1,1}) << std::endl;
      passedInternal = false;
    }
    if (hitMap.GetTotalCharge() != 3) {
      std::cout << " - FAILED " << std::endl << " | -> Expected total charge 3, got " << hitMap.GetTotalCharge() << std::endl; 
      passedInternal = false;
    }

    hitMap.FillCharge({0,0},100, simHitWrapper); // test adding charge to existing pixel

    if (hitMap.GetCharge({0,0}) != 101) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 101 at (0,0), got " << hitMap.GetCharge({0,0}) << std::endl;
      passedInternal = false;
    }
    if (hitMap.GetTotalCharge() != 103) {
      std::cout << " - FAILED " << std::endl << " | -> Expected total charge 103, got " << hitMap.GetTotalCharge() << std::endl; 
      passedInternal = false;
    }

    try {
      hitMap.FillCharge({-1,0},1, simHitWrapper);
      std::cout << " - FAILED " << std::endl << " | -> Expected out-of-bounds exception for pixel (-1,0)" << std::endl;
      passedInternal = false;
    }
    catch (const std::runtime_error& e) {}
    try {
      hitMap.FillCharge({0,-1},1, simHitWrapper);
      std::cout << " - FAILED " << std::endl << " | -> Expected out-of-bounds exception for pixel (0,-1)" << std::endl;
      passedInternal = false;
    }
    catch (const std::runtime_error& e) {}
    try {
      hitMap.FillCharge({10,0},1, simHitWrapper);
      std::cout << " - FAILED " << std::endl << " | -> Expected out-of-bounds exception for pixel (10,0)" << std::endl;
      passedInternal = false;
    }
    catch (const std::runtime_error& e) {}
    try {
      hitMap.FillCharge({0,10},1, simHitWrapper);
      std::cout << " - FAILED " << std::endl << " | -> Expected out-of-bounds exception for pixel (0,10)" << std::endl;
      passedInternal = false;
    }
    catch (const std::runtime_error& e) {}

    hitMap.Reset();

    if (hitMap.GetTotalCharge() != 0) {
      std::cout << " - FAILED " << std::endl << " | -> Expected total charge 0 after Reset(), got " << hitMap.GetTotalCharge() << std::endl; 
      passedInternal = false;
    }
    
    if (passedInternal)
      std::cout << " - PASSED" << std::endl;
    passed = passed && passedInternal;
  }

  if (passed)
    std::cout << " | All tests passed." << std::endl;
  return passed;
}

} // namespace VTXdigi_tools


