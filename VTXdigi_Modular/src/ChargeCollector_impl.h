// VTXdigi_Modular/src/ChargeCollector_impl.h
#pragma once

#include "TFile.h"
#include "TH1D.h"

#include "../include/VTXdigi_Modular.h"
#include <optional>

namespace VTXdigi_tools {

class SimHitWrapper; // forward-declaration for include/VTXdigi_tools.h
class HitMap;        // forward-declaration for include/VTXdigi_tools.h

constexpr double kMinPathCosTheta = 1e-6; // minimum cos(theta) below which a simHit is considered to be parallel to the
                                          // sensor surface. (theta = angle between simHit momentum and sensor normal)
constexpr float kPathLengthTolerance =
    1.001f; // tolerance factor for how much longer the computed path can be compared to the Geant4 path length.
// (If the computed path is longer than the Geant4 path, either the linear path approximation breaks down, or the
// particle begins or ends inside the sensor volume)

constexpr float kCCMapEntryMinimum =
    1.e-5f; // CCMap entries below this value are set to zero, to minimise unnecessary computations in hot loop.
// result is quite sensitive to this, so choose carefully. 1e-5 seems to be a good compromise between accuracy and
// performance for the TPSCo 65nm CIS charge collection maps

/** @brief Computes & then holds position & information about a simHits path through the sensor */
struct Path {
  dd4hep::rec::Vector3D entry;  // in local coordinates, in mm
  dd4hep::rec::Vector3D travel; // in local coordinates, in mm
  float lengthG4;
  float charge; // simHit deposited charge, in e
  bool isValid = false;

  Path() = default;
  /** @brief Construct path information from a simHit and the sensor's transformation matrix
   * @note set Path::valid. false means the path would not intersect the sensor volume
   * @param w_offset the offset in the w-direction (local z-direction). Used to get the same effect as shifting the
   * CCMap the opposite direction. Charge shifted outside the sensor volume is silenty droppend. When in doubt use
   * w_offset = 0.0 */
  Path(const SimHitWrapper& simHit, const TGeoHMatrix& trafoMatrix, const VTXdigi_Modular& digitizer, double w_offset);

  std::vector<std::pair<float, dd4hep::rec::Vector3D>>
  SampleDepositions(TRandom3& randomGen, const float meanDepositionsPerUm, const TH1D& chargeSamplingHist) const;
};

/** @brief Compute the factors by which to clip a path along a given axis (to clip it to the sensor volume) */
std::array<double, 2> ComputePathClippingFactors(std::array<double, 2> t, const double entry_ax, const double travel_ax,
                                                 const double sensorLength_ax);

/* -- Charge collector algorithm: Charge collection map-based -- */

class ChargeCollectionMap {
  /* TODO: use sparse storage (instead of storing ALL CCMap entries, even though ~80% are empty (for TPSCo 65nm CIS)).
   * Make sure to keep each matrix contiguous, because we iterate over those. (maybe even move to iterator over matrix
   * instead of getting each weight individually). This will also improve performance in HitMap::FillCharge() because
   * ~80% less operations have to be performed. */

  std::array<int, 3> m_voxelCount;
  int m_matrixSize;
  int m_matrixSize_half;         // =(matrixSize-1)/2 (for better performance in hot loop)
  std::vector<float> m_matrices; // charge sharing matrices, one per in-pixel bin (flattened for better performance)
  double m_chargeCollectionDepthCenter = 0.0;
  const bool m_shiftCollectionDepthToCenter;

public:
  /** @brief Construct charge collection map from a file
   * @note Checks that parameters from the charge collection map file match those in the digitiser */
  ChargeCollectionMap(const std::string& ccmapFileName, const VTXdigi_Modular& digitizer); // load weights from file

  double GetChargeCollectionDepthCenter() const { return m_chargeCollectionDepthCenter; }

  /** @brief Set the charge sharing matrix for a specific in-pixel bin
   * @param weights A flat vector containing the matrix weights in row-major order (ie. row-by-row) (length must be
   * matrixSize*matrixSize) */
  void SetMatrix(const VoxelIndex& voxI, const std::vector<float>& weights);

  /** @brief Set all charge sharing matrices to the same values */
  void SetAllMatrices(const std::vector<float>& weights);

  /** @brief Access matrix as const reference */
  const std::vector<float>& GetMatrix(const VoxelIndex& voxI) const;

  /** @brief Access a specific weight in a charge sharing matrix */
  inline float GetWeight(const VoxelIndex& voxI, const int col, const int row) const {
    return m_matrices[FindIndex(voxI, col, row)];
  };

  inline int GetSize() const { return m_matrixSize; }
  inline int GetSizeHalf() const { return m_matrixSize_half; }

  inline int GetVoxelCount(int i) const { return m_voxelCount.at(i); }
  inline std::array<int, 3> GetVoxelCount() const { return m_voxelCount; }

  /** @brief Get the vertical shift to apply to the path in case the CCMap collects charge assyemetrically around the
  vertical center
  @note Returns a positive value to shift the path upwards (to +w). This is mathematically equivalent to shifting the
  CCMap downwards, which is needed if the charge collection centre is above w=0.
  @note Returns 0 if Gaudi property m_shiftCollectionDepthToCenter==false */
  inline double VerticalPathShift() const {
    if (m_shiftCollectionDepthToCenter)
      return m_chargeCollectionDepthCenter;
    return 0.0;
  }

private:
  /** @brief Convert 3D in-pixel bin indices and a matrix row/column to a flat index for m_matrices */
  int FindIndex(const VoxelIndex& voxI, const int col, const int row) const;

}; // class ChargeCollectionMap

class ChargeCollector_CCMap : public IChargeCollector {

  ChargeCollectionMap m_CCMap;
  std::unique_ptr<const TH1D> m_chargeSamplingHist; // histogram to sample deposition charges from.
  float m_meanDepositionsPerUm; // mean number of deposition clusters per um of path length in the sensor to sample
                                // deposition count from

public:
  explicit ChargeCollector_CCMap(const VTXdigi_Modular& digitizer);

  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix,
               TRandom3& randomGen) const override;

private:
  /** @brief Fill charge into the hitmap, according to a in-pixel bin's charge sharing weights */
  void DistributeVoxelCharge(HitMap& hitMap, const float charge, const PixelIndex& pixI, const VoxelIndex& voxI,
                             const SimHitWrapper& simHit) const;
};

/* -- Charge collector algorithms for debugging -- */

class ChargeCollector_SinglePixel : public IChargeCollector {
public:
  explicit ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix,
               TRandom3& randomGen) const override;
};

class ChargeCollector_Debug : public IChargeCollector {
public:
  explicit ChargeCollector_Debug(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix,
               TRandom3& randomGen) const override;
};

/* -- Charge collector algorithm: Propagation-based -- */

class ChargeCollector_Drift : public IChargeCollector {
public:
  explicit ChargeCollector_Drift(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix,
               TRandom3& randomGen) const override;
};

} // namespace VTXdigi_tools
