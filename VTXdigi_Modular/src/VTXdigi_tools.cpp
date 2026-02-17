#include "VTXdigi_tools.h"
#include <iostream>

namespace VTXdigi_tools {


SimHitWrapper::SimHitWrapper(edm4hep::SimTrackerHit simHit, const dd4hep::rec::SurfaceMap* surfaceMap, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) : m_simHit(simHit) {
  m_cellID = GetCellID_short(simHit);
  const auto surfaceIt = surfaceMap->find(m_cellID);
  if (surfaceIt == surfaceMap->end())
    throw std::runtime_error("VTXdigi_Allpix2::HitInfo constructor: Could not find SimSurface for this hit's (reduced) cellID: " + std::to_string(m_cellID));
  // m_surface = std::make_shared<dd4hep::rec::ISurface>(surfaceIt->second);
  m_surface = surfaceIt->second;

  const float chargePerkeV = 273.97f; // in electrons, for silicon (1 eh-pair ~ 3.65 eV)
  m_charge = static_cast<int>(simHit.getEDep() * (dd4hep::GeV / dd4hep::keV) * chargePerkeV); // convert energy deposit (in keV) to number of electrons 

  m_layerNumber = GetLayer(m_cellID, cellIdDecoder);
} // Hit::Hit()

void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept {
  if (&a == &b) 
    return;

  using std::swap;
  swap(a.m_simHit, b.m_simHit);
  swap(a.m_surface, b.m_surface);
  swap(a.m_cellID, b.m_cellID);
  swap(a.m_charge, b.m_charge);
  swap(a.m_layerNumber, b.m_layerNumber);
} // swap(Hit&, Hit&)

/* -- helpers -- */

void CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitLinks, const dd4hep::rec::Vector3D& position, const int charge) {
  const float chargePerkeV = 273.97f; // in electrons, for silicon (1 eh-pair ~ 3.65 eV)

  auto digiHit = digiHits.create();

  digiHit.setCellID(simHit.getCellID());
  digiHit.setEDep(charge / chargePerkeV); // convert e- to keV
  digiHit.setPosition(ConvertVector(position));
  digiHit.setTime(simHit.getTime());
  
  auto digiHitLink = digiHitLinks.create();
  digiHitLink.setFrom(digiHit);
  digiHitLink.setTo(simHit);
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

dd4hep::DDSegmentation::CellID GetCellID_short(const edm4hep::SimTrackerHit& simHit) {
  /* Mask removes segmentation bits, now cellID is unique for each sensor */
  std::uint64_t m_mask = (static_cast<std::uint64_t>(1) << 32) - 1;
  return simHit.getCellID() & m_mask;
}

int GetLayer(const dd4hep::DDSegmentation::CellID& cellID, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) {
  return static_cast<int>(cellIdDecoder->get(cellID, "layer"));
} 
int GetLayer(const edm4hep::SimTrackerHit& simHit, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder) {
  return GetLayer(GetCellID_short(simHit), cellIdDecoder);
}

/* -- Binning things -- */

int ComputeBinIndex(float x, float binX0, float binWidth, int binN) {
  /** Get the bin index for a given x value
   *  binX0 is the lower edge of the first bin
   *  binWidth is the width of the bins
   *  binN is the number of bins
   *  return -1 if x is out of range
   */

  if (binN <= 0) throw GaudiException("ComputeBinIndex: binN must be positive", "VTXdigi_Allpix2::ComputeBinIndex()", StatusCode::FAILURE);
  if (binWidth <= 0.0) throw GaudiException("ComputeBinIndex: binWidth must be positive", "VTXdigi_Allpix2::ComputeBinIndex()", StatusCode::FAILURE);

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

std::array<int, 3> ComputeInPixelIndices(const dd4hep::rec::Vector3D& pos, const std::array<size_t, 3> binCount, const std::pair<float, float> pixelPitch, const std::pair<float, float> sensorLength, const float sensorThickness) {
  int j_u, j_v, j_w;

  const float posShifted_u = pos.x() + 0.5 * sensorLength.first; // shift to [0, length_u]
  if (posShifted_u < 0.0 || posShifted_u > sensorLength.first) {
    j_u = -1; // out of bounds
  }
  else {
    float posInPixel_u = std::fmod(posShifted_u,  pixelPitch.first);
    if (posInPixel_u < 0.0) posInPixel_u +=  pixelPitch.first; // ensure positive remainder
    j_u = ComputeBinIndex(posInPixel_u, 0.0,  pixelPitch.first / binCount[0], binCount[0]);
  }

  const float posShifted_v = pos.y() + 0.5 * sensorLength.second;
  if (posShifted_v < 0.0 || posShifted_v > sensorLength.second) {
    j_v = -1; // out of bounds
  }
  else {
    float posInPixel_v = std::fmod(posShifted_v, pixelPitch.second);
    if (posInPixel_v < 0.0) posInPixel_v += pixelPitch.second;
    j_v = ComputeBinIndex(posInPixel_v, 0.0, pixelPitch.second / binCount[1], binCount[1]);
  }

  // vertical (w) binning: shift to [0, thickness]
  const float posShifted_w = pos.z() + 0.5 * sensorThickness;
  j_w = ComputeBinIndex(posShifted_w, 0.0, sensorThickness / binCount[2], binCount[2]); // no fmod, so out-of-bounds is caught

  return {j_u, j_v, j_w};
} // ComputeInPixelIndices()

dd4hep::rec::Vector3D ComputePixelPos_local(const std::pair<int, int> pixelIndex, const std::pair<float, float> sensorLength,  const std::pair<float, float> pixelPitch, float depletedRegionDepthCenter) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */
  
  float u = (static_cast<float>(pixelIndex.first) + 0.5f) * pixelPitch.first - 0.5f * sensorLength.first; // in mm
  float v = (static_cast<float>(pixelIndex.second) + 0.5f) * pixelPitch.second - 0.5f * sensorLength.second;
  float w = depletedRegionDepthCenter;
  
  return dd4hep::rec::Vector3D(u, v, w); 
}

dd4hep::rec::Vector3D ComputePixelPos_local(const std::pair<int, int> pixelIndex, const std::pair<float, float> sensorLength, const std::pair<float, float> pixelPitch) {
  return ComputePixelPos_local(pixelIndex, sensorLength, pixelPitch, 0.f);
}

dd4hep::rec::Vector3D ComputePos_local(const std::pair<float, float> index, const std::pair<float, float> sensorLength,  const std::pair<float, float> pixelPitch, float depletedRegionDepthCenter) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */
  
  float u = (index.first + 0.5f) * pixelPitch.first - 0.5f * sensorLength.first; // in mm. Add 0.5*pixelPitch to shift from pixel edge to center, since index 0 is defined as the center of the pixel.
  float v = (index.second + 0.5f) * pixelPitch.second - 0.5f * sensorLength.second;
  float w = depletedRegionDepthCenter;
  
  return dd4hep::rec::Vector3D(u, v, w); 
}

dd4hep::rec::Vector3D ComputePos_local(const std::pair<float, float> index, const std::pair<float, float> sensorLength,  const std::pair<float, float> pixelPitch) {
  return ComputePos_local(index, sensorLength, pixelPitch, 0.f);
}


/* -- HitMap -- */

inline uint64_t PixToKey(std::pair<int, int> pix) {
  return (static_cast<uint64_t>(pix.first) << 32) | static_cast<uint32_t>(pix.second);
}
inline std::pair<int, int> PixFromKey(uint64_t key) {
  const int i_u = static_cast<int>(key >> 32);
  const int i_v = static_cast<int>(key & 0xFFFFFFFF);
  return std::make_pair(i_u, i_v);
}

HitMap::HitMap(std::pair<size_t, size_t> pixelCount) : m_pixCount(pixelCount) {
  const int inverseOccupancy = 2000; // assume occupancy, 5e-4 is quite conservative for Z-run
  m_pixels.reserve(pixelCount.first * pixelCount.second / inverseOccupancy); // avoid too many reallocations
}

void HitMap::FillCharge(std::pair<int, int> pix, int charge) {
  FillCharge(pix, charge, nullptr);
}

void HitMap::FillCharge(std::pair<int, int> pix, int charge, const edm4hep::SimTrackerHit* simHit) {
  if (_OutOfBounds(pix)) {
    throw std::runtime_error("HitMap::FillCharge: pixel i_u or i_v ( " + std::to_string(pix.first) + ", " + std::to_string(pix.second) + ") out of range");
  }
  const u_int64_t key = PixToKey(pix);
  m_pixels.try_emplace(key, pix); // does nothing if pixel already exists, otherwise creates it with default charge 0
  m_pixels[key].charge += charge;
  m_pixels[key].simHits.push_back(simHit); 
}

int HitMap::GetCharge(std::pair<int, int> pix) const {
  if (_OutOfBounds(pix)) {
    throw std::runtime_error("HitMap::GetCharge: pixel i_u or i_v ( " + std::to_string(pix.first) + ", " + std::to_string(pix.second) + ") out of range");
  }
  auto it = m_pixels.find(PixToKey(pix));
  if (it != m_pixels.end()) {
    return it->second.charge;
  }
  return 0; // if pixel not found, charge is 0
}

int HitMap::GetTotalCharge() const {
  int totalCharge = 0;
  for (const auto& [index, pixHit] : m_pixels) {
    totalCharge += pixHit.charge;
  }
  return totalCharge;
}

std::unordered_map<uint64_t, Pixel>& HitMap::Hits() {
  return m_pixels;
}

inline int HitMap::GetTotalPixelsWithCharge() const {
  return m_pixels.size();
}

inline void HitMap::Reset() {
  m_pixels.clear();
}

inline bool HitMap::_OutOfBounds(std::pair<int, int> pix) const { 
  return (
    pix.first < 0
    || pix.first >= static_cast<int>(m_pixCount.first)
    || pix.second < 0
    || pix.second >= static_cast<int>(m_pixCount.second)
  );
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

  std::cout << " | VTXdigi_tools::ComputePos_local()";
  {
    bool passedInternal = true;

    const std::pair<float, float> pixelPitch = { 1.0, 2.0 };
    const std::pair<float, float> sensorLength = { 10.0, 20.0 }; // 10 x 10 pixels

    std::array<std::pair<float, float>, 5> inputs = {{{0., 0.}}};
    inputs = {{
      {0., 0.},
      {4.5, 4.5},
      {9.0, 9.0},
      {-0.5, -0.5},
      {9.5, 9.5},
    }};
    std::array<dd4hep::rec::Vector3D, inputs.size()> expectedOutputs = {{dd4hep::rec::Vector3D(0.f, 0.f, 0.f)}};
    expectedOutputs = {{
      dd4hep::rec::Vector3D( -4.5, -9.0, 0.0 ),
      dd4hep::rec::Vector3D( 0., 0., 0. ),
      dd4hep::rec::Vector3D( 4.5, 9.0, 0.0 ),
      dd4hep::rec::Vector3D( -5.0, -10.0, 0.0 ),
      dd4hep::rec::Vector3D( 5.0, 10.0, 0.0 ),
    }};
    
    for (size_t i = 0; i < inputs.size(); ++i) {
      dd4hep::rec::Vector3D result = ComputePos_local(inputs.at(i), sensorLength, pixelPitch);
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

    const std::array<size_t, 3> binCount = { 4, 4, 4 };
    const std::pair<float, float> pixelPitch = { 1.0f, 2.0f }; // bin-width is 0.25 
    const std::pair<float, float> sensorLength = { 10.0f, 20.0f };
    const float sensorThickness = 1.f;

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
      std::array<int, 3> result = ComputeInPixelIndices(inputs.at(i), binCount, pixelPitch, sensorLength, sensorThickness);
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

    HitMap hitMap({10,10});

    hitMap.FillCharge({0,0},1);
    hitMap.FillCharge({1,1},2);

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

    hitMap.FillCharge({0,0},100); // test adding charge to existing pixel

    if (hitMap.GetCharge({0,0}) != 101) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 101 at (0,0), got " << hitMap.GetCharge({0,0}) << std::endl;
      passedInternal = false;
    }
    if (hitMap.GetTotalCharge() != 103) {
      std::cout << " - FAILED " << std::endl << " | -> Expected total charge 103, got " << hitMap.GetTotalCharge() << std::endl; 
      passedInternal = false;
    }

    try {
      hitMap.FillCharge({-1,0},1);
      std::cout << " - FAILED " << std::endl << " | -> Expected out-of-bounds exception for pixel (-1,0)" << std::endl;
      passedInternal = false;
    }
    catch (const std::runtime_error& e) {}
    try {
      hitMap.FillCharge({0,-1},1);
      std::cout << " - FAILED " << std::endl << " | -> Expected out-of-bounds exception for pixel (0,-1)" << std::endl;
      passedInternal = false;
    }
    catch (const std::runtime_error& e) {}
    try {
      hitMap.FillCharge({10,0},1);
      std::cout << " - FAILED " << std::endl << " | -> Expected out-of-bounds exception for pixel (10,0)" << std::endl;
      passedInternal = false;
    }
    catch (const std::runtime_error& e) {}
    try {
      hitMap.FillCharge({0,10},1);
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

  // std::cout << " | VTXdigi_tools::ComputePixelCenter_Local()";
  // {
  //   bool passedInternal = true;

  //   const double thickness_inner=0.5, thickness_outer=0.5;
  //   const dd4hep::rec::Vector3D u_val(1,0,0), v_val(0,1,0), n_val(0,0,1), origin_val(0,0,0);
  //   dd4hep::Box box(0.5,1.0,0.05); // half-lengths in cm
  //   dd4hep::Material air = dd4hep::Material("Air");
  //   dd4hep::Volume volume = dd4hep::Volume("dummyVol", box, air);
  //   const int id_val = 0;
  //   std::unique_ptr<dd4hep::rec::ISurface> surface = std::make_unique<dd4hep::rec::VolPlaneImpl>(dd4hep::rec::SurfaceType::Plane, thickness_inner, thickness_outer, u_val, v_val , n_val, origin_val, volume, id_val );    

  //   std::cout << "length along u: " << surface->length_along_u() << " mm, length along v: " << surface->length_along_v() << " mm" << std::endl;

  //   const std::pair<float, float> pixelPitch = {1.0f, 2.0f};

  //   std::pair<int, int> pixelIndex = {0,0};
  //   dd4hep::rec::Vector3D result = ComputePixelCenter_Local(pixelIndex, *surface, pixelPitch, 0.f);
  //   dd4hep::rec::Vector3D expected( -4.5f, -9.0f, 0.f );
  //   if (result != expected) {
  //     std::cout << " - FAILED " << std::endl << " | -> Expected pixel center at (" << expected.x() << ", " << expected.y() << ", " << expected.z() << "), got (" << result.x() << ", " << result.y() << ", " << result.z() << ")" << std::endl;
  //     passedInternal = false;
  //   }

    // delete surface;

  //   if (passedInternal)
  //     std::cout << " - PASSED" << std::endl;
  //   passed = passed && passedInternal;
  // }

  if (passed)
    std::cout << " | All tests passed." << std::endl;
  return passed;
}

} // namespace VTXdigi_tools


