#pragma once

#include <vector>
#include <numeric>
#include <tuple>

#include "Gaudi/Property.h"
#include "GaudiKernel/RndmGenerators.h"

#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "DDRec/Material.h"

#include "DD4hep/Objects.h"
#include "DD4hep/Volumes.h"
#include "DD4hep/Shapes.h"
#include "DD4hep/DetElement.h"

#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include <queue>

namespace VTXdigi_tools {

/* -- Tool tests -- */

bool ToolTest();


/* -- SimHitWrapper -- */


/** @brief Class to contain all information about a simTrackerHit that is needed for the digitization.
 * @note this is where the simTrackerHit is actually stored, everything else (pixelHit / cluster) will store pointers to this. */
class SimHitWrapper {
  edm4hep::SimTrackerHit m_simTrackerHit;
  dd4hep::DDSegmentation::CellID m_cellID; // without segmentation bits
  float m_charge;
  int m_layerNumber;
  mutable dd4hep::rec::Vector3D m_truthPos; // truth position of the simHit, in local coordinates. Can be adjusted by the charge collection algorithm. Only used for histogramming after filling the hits

public:
  SimHitWrapper(edm4hep::SimTrackerHit simTrackerHit, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder);
  SimHitWrapper(const SimHitWrapper& other) = default;
  SimHitWrapper(SimHitWrapper&& other) = default;
  SimHitWrapper() = default;

  inline void SetTruthPos(const dd4hep::rec::Vector3D& pos) const { m_truthPos = pos; } // only used for histogramming after filling the hits, so not really a problem that this is mutable

  friend void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept;
  inline const edm4hep::SimTrackerHit* hitPtr() const { return &m_simTrackerHit; }
  /** @brief Access the truth position of the simHit in local coordinates
   * @note might have been adjusted by ChargeCollector::FillHit() to account for charge collection effects. */
  inline const dd4hep::rec::Vector3D truthPos() const { return m_truthPos; }

  inline dd4hep::DDSegmentation::CellID cellID() const { return m_cellID; }
  inline float charge() const { return m_charge; }
  inline int layer() const { return m_layerNumber; }
};

void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept;

/* -- HitMap (and everything we need for it to work) -- */

/** @brief A pixel in the hit map. Can have multiple contributing simHits */
struct Pixel {
  float charge;
  std::unordered_set<const SimHitWrapper*> simHits;
  std::pair<int, int> index; // This info is saved in (a) the map key, and (b) here inside the Pixel object. This is inefficient. But it makes the code a bit nicer not having to pass the index around separately.

  Pixel(std::pair<int, int> pix) : charge(0.f), index(pix) {
    simHits.reserve(2); // avoid too many reallocations, will rarely see more than 2 simTrackerHits contributing to the same pixel
  }
  Pixel() : charge(0.f), index({-1, -1}) {
    simHits.reserve(2); // avoid too many reallocations, will rarely see more than 2 simTrackerHits contributing to the same pixel
  }
};

struct Hash_PairInt {
  size_t operator()(const std::pair<int,int>& i_uv) const noexcept {
    return (static_cast<uint64_t>(i_uv.first) << 32) ^ static_cast<uint32_t>(i_uv.second);
  }
};

using PixelMap = std::unordered_map<std::pair<int, int>, Pixel, Hash_PairInt>;

/** @brief HitMap of all pixel hits on a sensor
 * @note uses a std::unordered_map to only store pixels that have charge, which is more memory efficient for large pixel counts and low occupancy. */
class HitMap {
  /* I tried implementing a vector that contains every pixel, but for large pixel counts this is very memory-inefficient. Instead, this class uses a std::map to only store pixels that have charge. This is more memory efficient for sparse hits, with O(1) simHit/sensor/event this is at least a factor 100 faster that the vector approach. Maybe not the case to ttbar run with O(0.1%) pixel occupancy. */
  
  PixelMap m_pixels; // hit pixels, stored by value
  std::pair<size_t, size_t> m_pixCount; // size of the sensor in pixels

public:
  HitMap(std::pair<size_t, size_t> pixCount);

  /** @brief Add charge and a simHit to a pixel */
  void FillCharge(std::pair<int, int> i_uv, float charge, const SimHitWrapper& simHitWrapper);

  void ApplyChargeSmearing(const Rndm::Numbers& rndm_charge);
  void ApplyThreshold(const float threshold);

  /** @brief Get one pixel's collected charge */
  float GetCharge(std::pair<int, int> i_uv) const;

  /** @brief Get the total charge across all pixels */
  float GetTotalCharge() const;

  /** @brief Return a const reference to the internal pixel map */
  inline const PixelMap& Hits() const { return m_pixels; };

  /** @brief Return the number of pixels with charge */
  inline int GetTotalPixelsWithCharge() const { return m_pixels.size(); };

  inline void Reset() { m_pixels.clear(); };

private:
  /** @brief Returns true if the pixel is out of bounds */
  inline bool _OutOfBounds(std::pair<int, int> i_uv) const;
}; // class HitMap

/* -- Clusterization -- */

/** @brief Contains a clusters charge, pointers to all pixels in the cluster, and pointers to all contributing simHits */
struct Cluster {
  std::vector<const Pixel*> pixels; // pointer to pixels in this cluster (stored by value in HitMap::m_pixels)
  std::unordered_set<const SimHitWrapper*> simHits;
  float charge = 0.f;
};

/** @brief Compute the center position of a cluster, computed via charge-weighed center of gravity */
std::pair<float, float> ComputeClusterPos_Weighted(const Cluster& cluster);

/** @brief Get the indices of all direct neighbors of a pixel */
std::array<std::pair<int, int>, 4> GetDirectNeighbors(const std::pair<int, int>& i_uv);

/** @brief Clusterize pixels in a HitMap, using direct neighbors */
std::vector<Cluster> Clusterize_NextNeighbors(const HitMap& hitMap);

/** @brief Clusterize pixels in a HitMap, without any clustering (ie. each pixel is its own cluster) */
std::vector<Cluster> Clusterize_NoClustering(const HitMap& hitMap);

/* -- helpers -- */

/** @brief Convert a dd4hep::rec::Vector3D to edm4hep::Vector3d and vice-versa */
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec);
/** @copydoc dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec) */
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3f vec);
/** @copydoc dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec) */
edm4hep::Vector3d ConvertVector(dd4hep::rec::Vector3D vec);

/** @brief Compute the transformation matrix from global detector to local sensor frame for a given sensor cellID */
TGeoHMatrix ComputeSensorTrafoMatrix(const dd4hep::DDSegmentation::CellID& cellID, const dd4hep::VolumeManager& volumeManager, const TGeoRotation& sensorNormalRotation);

/** @brief Transform a position from global detector coordinates to local sensor coordinates, using the sensor transformation matrix */
dd4hep::rec::Vector3D GlobalToLocal(const dd4hep::rec::Vector3D& global, const TGeoHMatrix& M);
/** @brief Transform a position from local sensor coordinates to global detector coordinates, using the sensor transformation matrix */
dd4hep::rec::Vector3D LocalToGlobal(const dd4hep::rec::Vector3D& local, const TGeoHMatrix& M);

dd4hep::DDSegmentation::CellID GetCellID_short(const edm4hep::SimTrackerHit& simTrackerHit);
dd4hep::DDSegmentation::CellID GetCellID_short(const dd4hep::DDSegmentation::CellID& cellID);

int GetLayer(const dd4hep::DDSegmentation::CellID& cellID, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder);
int GetLayer(const edm4hep::SimTrackerHit& simTrackerHit, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder);

/* -- Binning tools -- */

/** @brief Given a histogram definition (x0, binWidth, nBins) and a value x, compute the bin index i in which x falls.
 * @return Int, -1 if x is out of range.
 * @note Bins are 0-indexed (vs ROOT's 1-indexing) */
int ComputeBinIndex(float x, float binX0, float binWidth, int binN);

/** @brief Compute the pixel indices (i_u, i_v) for a given (local) position inside the sensor */
std::pair<int, int> ComputePixelIndices(const dd4hep::rec::Vector3D& pos, const std::pair<float, float> pixelPitch, const std::pair<size_t, size_t> pixelCount);

/** @brief Compute the in-pixel indices (j_u, j_v, j_w) for a given (local) position inside the pixel and layer index
 *  @note Assumption: each layer has only 1 type if sensor */
std::array<int, 3> ComputeInPixelIndices(const dd4hep::rec::Vector3D& pos, const std::array<int, 3>& binCount, const std::pair<float, float>& pixelPitch, const std::array<float, 3>& sensorDimensions);

/** @brief Compute the center position of a given pixel (i_u,i_v) in sensor-local coordinates (u,v,w) 
 * 
 * @note The w coordinate is set to depletedRegionDepthCenter. 0 for center, +25 for sensor surface, +20 for TPSCo 65nm maps.
*/
dd4hep::rec::Vector3D ComputePosFromPixIndex_local(const std::pair<int, int> pixelIndex, const std::pair<float, float> sensorLength,  const std::pair<float, float> pixelPitch, float depletedRegionDepthCenter);
/** @brief Compute the center position of a given pixel (i_u,i_v) in sensor-local coordinates (u,v,0) */
dd4hep::rec::Vector3D ComputePosFromPixIndex_local(const std::pair<int, int> pixelIndex, const std::pair<float, float> sensorLength, const std::pair<float, float> pixelPitch);

/** @brief Compute the position of a given (pixel-)index (i_u,i_v) in sensor-local coordinates (u,v,0) .
 * @note index 0 indicates the center of the pixel, index -0.5 the lower edge and +0.5 the upper edge.  
 * @note Does not check if the position is within the sensor bounds!
 * @note The w coordinate is set to depletedRegionDepthCenter. 0 for center, +25 for sensor surface, +20 for TPSCo 65nm maps.
 */
dd4hep::rec::Vector3D ComputePosFromPixIndex_local(const std::pair<float, float> index, const std::pair<float, float> sensorLength,  const std::pair<float, float> pixelPitch, float depletedRegionDepthCenter);
/** @brief Compute the position of a given (pixel-)index (i_u,i_v) in sensor-local coordinates (u,v,0) .
 * @note index 0 indicates the center of the pixel, index -0.5 the lower edge and +0.5 the upper edge.  
 * @note Does not check if the position is within the sensor bounds!
 */
dd4hep::rec::Vector3D ComputePosFromPixIndex_local(const std::pair<float, float> index, const std::pair<float, float> sensorLength, const std::pair<float, float> pixelPitch);

} // namespace VTXdigi_tools