#pragma once

#include <memory>
#include "../include/VTXdigi_Modular.h"

namespace VTXdigi_tools {

class SimHitWrapper; // forward-declare things in include/VTXdigi_tools.h
class HitMap;

using ::VTXdigi_Modular;

/** @brief Holds position & information about path through the sensor */
struct Path {
  Path(const std::shared_ptr<edm4hep::SimTrackerHit> simTrackerHit, const TGeoHMatrix& trafoMatrix, const VTXdigi_Modular& digitizer);
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

using Index_pix = std::pair<int, int>;
using Index_inPix = std::array<int, 3>;

class LookupTable {
  Index_inPix m_binCount;
  int m_matrixSize;
  std::vector<std::vector<float>> m_matrices; // charge sharing matrices, one per in-pixel bin
  /* TODO: improve performance by completely flattening the m_matrices vector? Propably not the bottleneck though... */

public:

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
   * @param i_uv {i_u, i_v} In-pixel indices
   * @param j_uvw {j_u, j_v, j_w} Column and row of the matrix entry to access (go from -(matrixSize-1)/2 to +(matrixSize-1)/2) */
  float GetWeight(const Index_inPix& j_uvw, const int col, const int row) const;

  inline int GetSize() const { return m_matrixSize; }
  inline int GetBinCountU() const { return m_binCount.at(0); }
  inline int GetBinCountV() const { return m_binCount.at(1); }
  inline int GetBinCountW() const { return m_binCount.at(2); }

private:

  int _FindIndex (const Index_inPix& j_uvw) const;
}; // class LookupTable


class ChargeCollector_LUT : public IChargeCollector {

  LookupTable m_LUT;

public:
  explicit ChargeCollector_LUT(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
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