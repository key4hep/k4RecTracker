#include "VTXdigi_tools.h"
#include <iostream>

namespace VTXdigi_tools {



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

std::array<int, 2> ComputePixelIndices(const dd4hep::rec::Vector3D& pos, const std::array<float, 2> pixelPitch, const std::array<size_t, 2> pixelCount) {
  
  const float length_u_half = 0.5 * pixelPitch.at(0) * pixelCount.at(0);
  int i_u = ComputeBinIndex(
    pos.x(),
    -length_u_half,
    pixelPitch.at(0),
    pixelCount.at(0));
    
  const float length_v_half = 0.5 * pixelPitch.at(1) * pixelCount.at(1);
  int i_v = ComputeBinIndex(
    pos.y(),
    -length_v_half,
    pixelPitch.at(1),
    pixelCount.at(1));

  return {i_u, i_v};
} // ComputePixelIndices()

std::array<int, 3> ComputeInPixelIndices(const dd4hep::rec::Vector3D& pos, const std::array<size_t, 3> binCount, const std::array<float, 2> pixelPitch, const std::array<float, 2> sensorLength, const float sensorThickness) {
  int j_u, j_v, j_w;

  const float posShifted_u = pos.x() + 0.5 * sensorLength.at(0); // shift to [0, length_u]
  if (posShifted_u < 0.0 || posShifted_u > sensorLength.at(0)) {
    j_u = -1; // out of bounds
  }
  else {
    float posInPixel_u = std::fmod(posShifted_u,  pixelPitch.at(0));
    if (posInPixel_u < 0.0) posInPixel_u +=  pixelPitch.at(0); // ensure positive remainder
    j_u = ComputeBinIndex(posInPixel_u, 0.0,  pixelPitch.at(0) / binCount[0], binCount[0]);
  }

  const float posShifted_v = pos.y() + 0.5 * sensorLength.at(1);
  if (posShifted_v < 0.0 || posShifted_v > sensorLength.at(1)) {
    j_v = -1; // out of bounds
  }
  else {
    float posInPixel_v = std::fmod(posShifted_v, pixelPitch.at(1));
    if (posInPixel_v < 0.0) posInPixel_v += pixelPitch.at(1);
    j_v = ComputeBinIndex(posInPixel_v, 0.0, pixelPitch.at(1) / binCount[1], binCount[1]);
  }

  // vertical (w) binning: shift to [0, thickness]
  const float posShifted_w = pos.z() + 0.5 * sensorThickness;
  j_w = ComputeBinIndex(posShifted_w, 0.0, sensorThickness / binCount[2], binCount[2]); // no fmod, so out-of-bounds is caught

  return {j_u, j_v, j_w};
} // ComputeInPixelIndices()



dd4hep::rec::Vector3D ComputePixelCenter_Local(const std::array<int, 2> pixelIndex, const dd4hep::rec::ISurface& surface,  const std::array<float, 2> pixelPitch, float depletedRegionDepthCenter) {
  /* returns the position of the center of pixel i_u, i_v in the local sensor frame */
  
  float length_u = surface.length_along_u() * 10; // convert to mm (works, checked 2025-10-17)
  float length_v = surface.length_along_v() * 10;
  
  float u = -0.5 * length_u + (pixelIndex[0] + 0.5) * pixelPitch[0]; // in mm
  float v = -0.5 * length_v + (pixelIndex[1] + 0.5) * pixelPitch[1];
  float w = depletedRegionDepthCenter;
  
  return dd4hep::rec::Vector3D(u, v, w); 
}

dd4hep::rec::Vector3D ComputePixelCenter_Local(const std::array<int, 2> pixelIndex, const dd4hep::rec::ISurface& surface, const std::array<float, 2> pixelPitch) {
  return ComputePixelCenter_Local(pixelIndex, surface, pixelPitch, 0.f);
}

/* -- PixelChargeMatrix -- */

PixelChargeMatrix::PixelChargeMatrix(int i_origin_u, int i_origin_v) : m_pixelCharge(m_initialSize * m_initialSize, 0.f), m_origin{ i_origin_u, i_origin_v } {
    // At first, the range is centered around the origin
  const int halfSize = (m_initialSize - 1) / 2;
  m_range_u[0] = m_origin[0] - halfSize;
  m_range_u[1] = m_origin[0] + halfSize;
  m_range_v[0] = m_origin[1] - halfSize;
  m_range_v[1] = m_origin[1] + halfSize;
}

void PixelChargeMatrix::Reset() {
  std::fill(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f);
}

inline std::tuple<int, int> PixelChargeMatrix::GetSize() const {
  return std::make_tuple(GetSize_u(), GetSize_v());
}
inline std::tuple<int, int> PixelChargeMatrix::GetOrigin() const {
  return std::make_tuple(m_origin[0], m_origin[1]);
}

float PixelChargeMatrix::GetTotalCharge() const {
  return std::accumulate(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f);
}

float PixelChargeMatrix::GetCharge(int i_u, int i_v) const {
  if (_OutOfBounds(i_u, i_v)) {
    throw std::runtime_error("PixelChargeMatrix::GetCharge: pixel i_u or i_v ( " + std::to_string(i_u) + ", " + std::to_string(i_v) + ") out of range");
  }
  return m_pixelCharge[_FindIndex(i_u, i_v)];
}

void PixelChargeMatrix::FillCharge(int i_u, int i_v, float charge) {
  if (_OutOfBounds(i_u, i_v)) {
    _ExpandMatrix(i_u, i_v);
  }
  m_pixelCharge[_FindIndex(i_u, i_v)] += charge;
}

inline int PixelChargeMatrix::_FindIndex(int i_u, int i_v) const {
  int i_u_rel = i_u - m_range_u[0];
  int i_v_rel = i_v - m_range_v[0];
  return i_u_rel + i_v_rel * GetSize_u();
}

inline bool PixelChargeMatrix::_OutOfBounds(int i_u, int i_v) const {
  return (
    i_u < m_range_u[0]
    || i_u > m_range_u[1]
    || i_v < m_range_v[0]
    || i_v > m_range_v[1]);
}

void PixelChargeMatrix::_ExpandMatrix(int i_u, int i_v) {
  /* Expand the matrix to include (i_u, i_v) and an excess of m_expansionStep pixels. */

  int rangeNew_u[2] = { m_range_u[0], m_range_u[1] }; 
  int rangeNew_v[2] = { m_range_v[0], m_range_v[1] };

  // TODO: optimise the expansion logic below, to expand to two directions at once if hit is close to corner ~ Jona 2025-10

  /* Expand it in each direction if necessary */
  if (i_u < m_range_u[0]) {
    rangeNew_u[0] = i_u - m_overExpansionStep;
  } else if (i_u > m_range_u[1]) {
    rangeNew_u[1] = i_u + m_overExpansionStep;
  }
  if (i_v < m_range_v[0]) {
    rangeNew_v[0] = i_v - m_overExpansionStep;
  } else if (i_v > m_range_v[1]) {
    rangeNew_v[1] = i_v + m_overExpansionStep;
  }

  const int newSizeU = rangeNew_u[1] - rangeNew_u[0] + 1;
  const int newSizeV = rangeNew_v[1] - rangeNew_v[0] + 1;

  std::vector<float> newPixelCharge(newSizeU * newSizeV, 0.f);

  /* Copy old charges into new vector, row by row */
  for (int row = 0; row < GetSize_v(); ++row) {
    std::copy(
      m_pixelCharge.begin() + row * GetSize_u(),
      m_pixelCharge.begin() + (row + 1) * GetSize_u(),
      newPixelCharge.begin() + (row + (m_range_v[0] - rangeNew_v[0])) * newSizeU + (m_range_u[0] - rangeNew_u[0])
    );
  }

  m_pixelCharge.swap(newPixelCharge);
  std::copy(rangeNew_u, rangeNew_u + 2, m_range_u);
  std::copy(rangeNew_v, rangeNew_v + 2, m_range_v);
} // _ExpandMatrix()


/* -- SensorChargeMatrix -- */

SensorChargeMatrix::SensorChargeMatrix(std::array<size_t, 2> pixelCount) : m_sensorCharge(pixelCount.at(0) * pixelCount.at(1), 0), m_pixelCount(pixelCount) {}

void SensorChargeMatrix::FillCharge(std::array<int, 2> pixel, int charge)
{
  if (_OutOfBounds(pixel)) {
    throw std::runtime_error("SensorChargeMatrix::FillCharge: pixel i_u or i_v ( " + std::to_string(pixel[0]) + ", " + std::to_string(pixel[1]) + ") out of range");
  }
  m_sensorCharge[_FindIndex(pixel)] += charge;
}

int SensorChargeMatrix::GetTotalCharge() const {
  return std::accumulate(m_sensorCharge.begin(), m_sensorCharge.end(), 0);
}

int SensorChargeMatrix::GetCharge(std::array<int, 2> pixel) const {
  if (_OutOfBounds(pixel)) {
    throw std::runtime_error("SensorChargeMatrix::GetCharge: pixel i_u or i_v ( " + std::to_string(pixel[0]) + ", " + std::to_string(pixel[1]) + ") out of range");
  }
  return m_sensorCharge[_FindIndex(pixel)];
}

void SensorChargeMatrix::Reset() {
  std::fill(m_sensorCharge.begin(), m_sensorCharge.end(), 0);
}

inline bool SensorChargeMatrix::_OutOfBounds(std::array<int, 2> pixel) const {
  return (
    pixel[0] < 0
    || pixel[0] >= static_cast<int>(m_pixelCount.at(0))
    || pixel[1] < 0
    || pixel[1] >= static_cast<int>(m_pixelCount.at(1))
  );
}

inline int SensorChargeMatrix::_FindIndex(std::array<int, 2> pixel) const {
  return pixel[0] + pixel[1] * m_pixelCount.at(0);
}


/* Tool tests */

bool ToolTest() {
  std::cout << " | Running VTXdigi tool tests" << std::endl;
  
  std::cout << " | VTXdigi_tools::ComputeBinIndex()";
  bool passed = true;
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

    const std::array<float, 2> pixelPitch = { 1.0f, 2.0f };
    const std::array<size_t, 2> pixelCount = { 10, 10 };

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
    std::array<std::array<int, 2>, inputs.size()> expectedOutputs = {{{0, 0}}};
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
      std::array<int, 2> result = ComputePixelIndices(inputs.at(i), pixelPitch, pixelCount);
      if (result != expectedOutputs[i]) {
        std::cout << " - FAILED " << std::endl << " | -> Expected pixel indices (" << expectedOutputs[i][0] << ", " << expectedOutputs[i][1] << ") for (u,v,w)=(" << inputs.at(i).x() << ", " << inputs.at(i).y() << ", " << inputs.at(i).z() << "), got (" << result[0] << ", " << result[1] << ")" << std::endl;
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
    const std::array<float, 2> pixelPitch = { 1.0f, 2.0f }; // bin-width is 0.25 
    const std::array<float, 2> sensorLength = { 10.0f, 20.0f };
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

  std::cout << " | VTXdigi_tools::PixelChargeMatrix()";
  {
    bool passedInternal = true;

    PixelChargeMatrix matrix(0,0);

    matrix.FillCharge(0,0,1.0f);
    matrix.FillCharge(1,1,2.0f);

    if (matrix.GetCharge(0,0) != 1.0f) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 1.0f at (0,0), got " << matrix.GetCharge(0,0) << std::endl;
      passedInternal = false;
    }
    if (matrix.GetCharge(1,1) != 2.0f) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 2.0f at (1,1), got " << matrix.GetCharge(1,1) << std::endl;
      passedInternal = false;
    }
    if (matrix.GetTotalCharge() != 3.0f) {
      std::cout << " - FAILED " << std::endl << " | -> Expected total charge 3.0f, got " << matrix.GetTotalCharge() << std::endl; 
      passedInternal = false;
    }

    matrix.Reset();
    if (matrix.GetTotalCharge() != 0.0f) {
      std::cout << " - FAILED " << std::endl << " | -> Expected total charge 0.0f after Reset(), got " << matrix.GetTotalCharge() << std::endl; 
      passedInternal = false;
    }

    std::tuple<int, int> sizeBefore = matrix.GetSize();
    matrix.FillCharge(20,20,1.5f); // should trigger expansion
    if (matrix.GetCharge(20,20) != 1.5f) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 1.5f at (20,20), got " << matrix.GetCharge(20,20) << std::endl;
      passedInternal = false;
    }
    std::tuple<int, int> sizeAfter = matrix.GetSize();
    if (std::get<0>(sizeAfter) <= std::get<0>(sizeBefore) || std::get<1>(sizeAfter) <= std::get<1>(sizeBefore)) {
      std::cout << " - FAILED " << std::endl << " | -> Expected matrix size to increase after expansion, before: (" << std::get<0>(sizeBefore) << ", " << std::get<1>(sizeBefore) << "), after: (" << std::get<0>(sizeAfter) << ", " << std::get<1>(sizeAfter) << ")" << std::endl;
      passedInternal = false;
    }
    
    if (passedInternal)
      std::cout << " - PASSED" << std::endl;
    passed = passed && passedInternal;
  }

  std::cout << " | VTXdigi_tools::SensorChargeMatrix()";
  {
    bool passedInternal = true;

    const std::array<size_t, 2> pixelCount = { 10, 10 };
    SensorChargeMatrix matrix(pixelCount);

    matrix.FillCharge({0,0}, 1);
    matrix.FillCharge({9,9}, 2);

    if (matrix.GetCharge({0,0}) != 1) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 1 at (0,0), got " << matrix.GetCharge({0,0}) << std::endl;
      passedInternal = false;
    }
    if (matrix.GetCharge({9,9}) != 2) {
      std::cout << " - FAILED " << std::endl << " | -> Expected charge 2 at (9,9), got " << matrix.GetCharge({9,9}) << std::endl;
      passedInternal = false;
    }
    if (matrix.GetTotalCharge() != 3) {
      std::cout << " - FAILED " << std::endl << " | -> Expected total charge 3, got " << matrix.GetTotalCharge() << std::endl; 
      passedInternal = false;
    }

    matrix.Reset();
    if (matrix.GetTotalCharge() != 0) {
      std::cout << " - FAILED " << std::endl << " | -> Expected total charge 0 after Reset(), got " << matrix.GetTotalCharge() << std::endl; 
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

  //   const std::array<float, 2> pixelPitch = {1.0f, 2.0f};

  //   std::array<int, 2> pixelIndex = {0,0};
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


