#pragma once

#include <memory>
#include "../include/VTXdigi_Modular.h"

namespace VTXdigi_tools {

class SimHitWrapper; // forward-declare things in include/VTXdigi_tools.h
class HitMap;

using ::VTXdigi_Modular;

using Index_pix = std::pair<int, int>;
using Index_inPix = std::array<int, 3>;

struct Index_segment {
  Index_pix i;
  Index_inPix j;

  inline bool operator==(const Index_segment& other) const {
    return (i == other.i) && (j == other.j);
  }
};

/** @brief Holds position & information about path through the sensor */
struct Path {

  /** @brief Construct path information from a simHit and the sensor's transformation matrix 
   * @note Only returns paths inside the sensor volume */
  Path(const edm4hep::SimTrackerHit* simTrackerHit, const TGeoHMatrix& trafoMatrix, const VTXdigi_Modular& digitizer);
  Path() = default;

  dd4hep::rec::Vector3D entry;
  dd4hep::rec::Vector3D travel;
  dd4hep::rec::Vector3D simPos;

  float length; // in mm
  float lengthG4;
  int nSegments;
  float charge;
};

std::pair<float, float> ComputePathClippingFactors(std::pair<float,float> t, const float entry_ax, const float travel_ax, const float sensorLength_ax, const VTXdigi_Modular& digitizer);



/* -- Charge collector algorithm: LUT-based -- */

struct VoxelHit { 
  Index_pix i;
  Index_inPix j;

  float t0, t1;
  float len;

  VoxelHit(Index_pix i_, Index_inPix j_, float len_) : i(i_), j(j_), len(len_) {}
};

inline float safe_div(float a, float b) {
  if (b == 0.f) {
    return std::numeric_limits<float>::infinity();
  }
  return a/b;
}


class LookupTable {
  Index_inPix m_binCount;
  int m_matrixSize;
  int m_matrixSize_half; // (matrixSize-1)/2
  // std::vector<std::vector<float>> m_matrices; // charge sharing matrices, one per in-pixel bin
  std::vector<float> m_matrices; // charge sharing matrices, one per in-pixel bin
  /* TODO: improve performance by completely flattening the m_matrices vector? Propably not the bottleneck though... */

public:

  /** @brief Construct lookup table from a file
   * @note Checks that parameters from the LUT file match those in the digitiser */
  LookupTable(const std::string& lutFileName, const VTXdigi_Modular& digitizer); // load weights from file
  LookupTable() = default;

  /** @brief Set the charge sharing matrix for a specific in-pixel bin 
   * @param weights A flat vector containing the matrix weights in row-major order (ie. row-by-row) (length must be matrixSize*matrixSize)
   */
  void SetMatrix(const Index_inPix& j_uvw, const std::vector<float>& weights);

  /** @brief Set all charge sharing matrices to the same values  */
  void SetAllMatrices(const std::vector<float>& weights);

  /** @brief Access matrix as const reference */
  const std::vector<float>& GetMatrix(const Index_inPix& j_uvw) const;

  /** @brief Access a specific weight in a charge sharing matrix
   * @note this seems to bottleneck the digitizer a lot */
  inline float GetWeight(const Index_inPix& j_uvw, const int col, const int row) const {
    return m_matrices[_FindIndex(j_uvw, col, row)];
  }; 

  inline int GetSize() const { return m_matrixSize; }
  inline int GetSizeHalf() const { return m_matrixSize_half; }

  inline int GetBinCount(int i) const { return m_binCount.at(i); }
  inline Index_inPix GetBinCount() const { return m_binCount; }
private:

  int _FindIndex (const Index_inPix& j, const int col, const int row) const;
}; // class LookupTable


class ChargeCollector_LUT : public IChargeCollector {

  LookupTable m_LUT;
  const float m_stepLength; // in mm

public:

  explicit ChargeCollector_LUT(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;

private:

  Index_segment ComputeSegmentIndices(const int step, const int stepCount, const Path& path) const;
  void DistributeSegmentCharge(HitMap& hitMap, const Index_segment& i_seg, const float charge, const int segmentsInBin, const edm4hep::SimTrackerHit* m_simTrackerHit) const;
};


/* -- Charge collector algorithms for debugging -- */

class ChargeCollector_SinglePixel : public IChargeCollector {
public:
  explicit ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};

class ChargeCollector_Debug : public IChargeCollector {
public:
  explicit ChargeCollector_Debug(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};



/* -- Charge collector algorithm: Propagation-based -- */

class ChargeCollector_Drift : public IChargeCollector {
public:
  explicit ChargeCollector_Drift(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};

} // namespace VTXdigi_tools