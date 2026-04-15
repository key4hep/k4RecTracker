// VTXdigi_Modular/src/ChargeCollector_impl.h
#pragma once

#include "../include/VTXdigi_Modular.h"

namespace VTXdigi_tools {
class SimHitWrapper; // forward-declaration for include/VTXdigi_tools.h
class HitMap; // forward-declaration for include/VTXdigi_tools.h

using Index_pix = std::pair<int, int>;
using Index_inPix = std::array<int, 3>;

constexpr float kPathLengthTolerance = 1.05f; // tolerance factor for how much longer the computed path can be compared to the Geant4 path length. If the computed path is longer than the Geant4 path, either the linear path approximation breaks down, or the particle begins or ends inside the sensor volume 
constexpr float kLutEntryMinimum = 1.e-3f; // LUT entries below this value are set to zero, to minimise unnecessary computations in hot loop. result is quite sensitive to this, so choose carefully. 1e-5 seems to be a good compromise between accuracy and performance for the TPSCo 65nm CIS LUT

/** @brief holds pixel indices i and in-pixel bin indices j (eg. for a segment) */
struct Index_segment {
  Index_pix i;
  Index_inPix j;

  inline bool operator==(const Index_segment& other) const {
    return (i == other.i) && (j == other.j);
  }
};

/** @brief Computes & then holds position & information about a simHits path through the sensor */
struct Path {
  Path() = default;

  dd4hep::rec::Vector3D entry; // in local coordinates, in mm
  dd4hep::rec::Vector3D travel; // in local coordinates, in mm
  dd4hep::rec::Vector3D simPos; // in local coordinates, in mm

  float length; // in mm
  float lengthG4;
  int nSegments;
  float charge;
};

/** @brief Compute the factors by which to clip a path along a given axis (to clip it to the sensor volume) */
std::pair<float, float> ComputePathClippingFactors(std::pair<float,float> t, const float entry_ax, const float travel_ax, const float sensorLength_ax);

/** @brief Construct path information from a simHit and the sensor's transformation matrix 
 * @note returns true if path is valid, false otherwise. false means the path would not intersect the sensor volume */
bool ConstructPath(Path& path, const SimHitWrapper& simHit, const TGeoHMatrix& trafoMatrix, const VTXdigi_Modular& digitizer);


/* -- Charge collector algorithm: LUT-based -- */

class LookupTable {
  /* TODO: use sparse storage (instead of storing ALL lut entries, even though ~80% are empty (for TPSCo 65nm CIS)). Make sure to keep each matrix contiguous, because we iterate over those. (maybe even move to iterator over matrix instead of getting each weight individually). This will also improve performance in HitMap::FillCharge() because ~80% less operations have to be performed. */

  Index_inPix m_binCount;
  int m_matrixSize;
  int m_matrixSize_half; // =(matrixSize-1)/2 (for better performance in hot loop)
  std::vector<float> m_matrices; // charge sharing matrices, one per in-pixel bin (flattened for better performance)

public:

  /** @brief Construct lookup table from a file
   * @note Checks that parameters from the LUT file match those in the digitiser */
  LookupTable(const std::string& lutFileName, const VTXdigi_Modular& digitizer); // load weights from file
  LookupTable() = default;

  /** @brief Set the charge sharing matrix for a specific in-pixel bin 
   * @param weights A flat vector containing the matrix weights in row-major order (ie. row-by-row) (length must be matrixSize*matrixSize) */
  void SetMatrix(const Index_inPix& j_uvw, const std::vector<float>& weights);

  /** @brief Set all charge sharing matrices to the same values */
  void SetAllMatrices(const std::vector<float>& weights);

  /** @brief Access matrix as const reference */
  const std::vector<float>& GetMatrix(const Index_inPix& j_uvw) const;

  /** @brief Access a specific weight in a charge sharing matrix */
  inline float GetWeight(const Index_inPix& j_uvw, const int col, const int row) const {
    return m_matrices[FindIndex(j_uvw, col, row)];
  }; 

  inline int GetSize() const { return m_matrixSize; }
  inline int GetSizeHalf() const { return m_matrixSize_half; }

  inline int GetBinCount(int i) const { return m_binCount.at(i); }
  inline Index_inPix GetBinCount() const { return m_binCount; }
private:

  /** @brief Convert 3D in-pixel bin indices and a matrix row/column to a flat index for m_matrices */
  int FindIndex (const Index_inPix& j, const int col, const int row) const;
}; // class LookupTable

class ChargeCollector_LUT : public IChargeCollector {

  LookupTable m_LUT;
  const float m_stepLength; // in mm

public:

  explicit ChargeCollector_LUT(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;

private:

  /** @brief Compute the pixel- and in-pixel bin indices for a given segment in a path */
  Index_segment ComputeSegmentIndices(const int step, const int stepCount, const Path& path) const;

  /** @brief Fill charge into the hitmap, according to a in-pixel bin's charge sharing weights */
  void DistributeSegmentCharge(HitMap& hitMap, const Index_segment& i_seg, const float charge, const int segmentsInBin, const SimHitWrapper& simHit) const;
  
  /** @brief Move the truth position in the simHitWrapper along the computed path to the depleted region center (as defined by the digitizer Gaudi property) 
   * @note This is necessary to get meaningful residual plots for sensors where charge is not collected across the whole sensor thickness. In TPSCo 65nm CIS the depletion region only extends ~10 um down into the sensor, so charges are only collected close to the upper sensor surface. For tracks passing through the sensor, the simHit position is defined at the sensor center. So for angled tracks, the residual (which is truth - reco) would be dominated by the track angle instead of the charge collection effects we want to study. By moving the truth position to the depleted region center, we get residuals that are dominated by charge collection effects.
   * @note Only used for plotting the residuals, does not change the output collections in any way.
   * @note If the path does not reach the depleted region center, the truth position is moved to the path end closest to the depleted region center */
  void MoveTruthPosition(const SimHitWrapper& simHit, const Path& path) const;
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