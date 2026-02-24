#pragma once

#include <vector>
#include <numeric>
#include <tuple>

#include "Gaudi/Property.h"

// #include "DDRec/ISurface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "DDRec/Material.h"

#include "DD4hep/Objects.h"
#include "DD4hep/Volumes.h"
#include "DD4hep/Shapes.h"
#include "DD4hep/DetElement.h"

// #include "edm4hep/SimTrackerHitCollection.h"
// #include "edm4hep/EventHeaderCollection.h" 
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include <queue>

namespace VTXdigi_tools {



/* -- Tool tests -- */

bool ToolTest();


/* -- SimHitWrapper -- */

class SimHitWrapper {
  /* any member that is added needs to be added to swap, too! */
  std::shared_ptr<edm4hep::SimTrackerHit> m_simTrackerHit;
  dd4hep::rec::ISurface* m_surface;
  
  dd4hep::DDSegmentation::CellID m_cellID; // cellID (without segmentation bits)
  float m_charge;
  int m_layerNumber;

public:
  SimHitWrapper(edm4hep::SimTrackerHit simTrackerHit, const dd4hep::rec::SurfaceMap* surfaceMap, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder);
  SimHitWrapper(const SimHitWrapper& other) = default;
  SimHitWrapper(SimHitWrapper&& other) = default;

  friend void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept;
  inline const std::shared_ptr<edm4hep::SimTrackerHit> hitPtr() const { return m_simTrackerHit; }
  inline dd4hep::rec::ISurface* surface() const { return m_surface; }

  inline dd4hep::DDSegmentation::CellID cellID() const { return m_cellID; }
  inline float charge() const { return m_charge; }
  inline int layer() const { return m_layerNumber; }
};

void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept;

/* -- HitMap (and everything we need for it to work) -- */

struct Pixel {
  float charge;
  std::unordered_set<std::shared_ptr<const edm4hep::SimTrackerHit>> simTrackerHits;
  std::pair<int, int> index; // This info is saved in (a) the map key, and (b) here inside the Pixel object. This is inefficient. But it makes the code a bit nicer not having to pass the index around separately.

  Pixel(std::pair<int, int> pix) : charge(0.f), index(pix) {
    simTrackerHits.reserve(2); // avoid too many reallocations, will rarely see more than 2 simTrackerHits contributing to the same pixel
  }
  Pixel() : charge(0.f), index({-1, -1}) {
    simTrackerHits.reserve(2); // avoid too many reallocations, will rarely see more than 2 simTrackerHits contributing to the same pixel
  }
};

struct PairHash {
  size_t operator()(const std::pair<int,int>& i_uv) const noexcept {
    return (static_cast<uint64_t>(i_uv.first) << 32) ^ static_cast<uint32_t>(i_uv.second);
  }
};

using PixelMap = std::unordered_map<std::pair<int, int>, std::shared_ptr<Pixel>, PairHash>;

class HitMap {
  /* I tried implementing a vector that contains every pixel, but for large pixel counts this is very memory-inefficient. Instead, this class uses a std::map to only store pixels that have charge. This is more memory efficient for sparse hits, with O(1) simHit/sensor/event this is at least a factor 100 faster that the vector approach. Maybe not the case to ttbar run with O(0.1%) pixel occupancy. */
  
  PixelMap m_pixels;
  std::pair<size_t, size_t> m_pixCount; // number of pixels in u and v direction

public:
  HitMap(std::pair<size_t, size_t> pixCount);

  void FillCharge(std::pair<int, int> i_uv, float charge, std::shared_ptr<const edm4hep::SimTrackerHit> simTrackerHit);
  float GetCharge(std::pair<int, int> i_uv) const;
  float GetTotalCharge() const;
  // inline const std::unordered_map<std::pair<int, int>, Pixel, PairHash>& Hits() const { return m_pixels; };
  // inline std::shared_ptr<const std::unordered_map<std::pair<int, int>, Pixel, PairHash>> HitsPtr() const { return std::make_shared<const std::unordered_map<std::pair<int, int>, Pixel, PairHash>>(m_pixels); };
  inline const PixelMap& Hits() const { return m_pixels; };
  inline int GetTotalPixelsWithCharge() const { return m_pixels.size(); };
  inline void Reset() { m_pixels.clear(); };

private:
  inline bool _OutOfBounds(std::pair<int, int> i_uv) const;
}; // class HitMap

/* -- Clusterization -- */

struct Cluster {
  std::vector<std::shared_ptr<Pixel>> pixels;
  std::unordered_set<std::shared_ptr<const edm4hep::SimTrackerHit>> simTrackerHits;
  float charge = 0.f;
};

std::pair<float, float> ComputeClusterPos_Weighted(const Cluster& cluster);

std::array<std::pair<int, int>, 4> GetDirectNeighbors(const std::pair<int, int>& i_uv);

std::vector<Cluster> Clusterize_NextNeighbors(const HitMap& hitMap);

std::vector<Cluster> Clusterize_NoClustering(const HitMap& hitMap);

/* -- helpers -- */

void CreateDigiHit(const edm4hep::SimTrackerHit& simTrackerHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitLinks, const dd4hep::rec::Vector3D& position, const float charge);

dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec);
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3f vec);
edm4hep::Vector3d ConvertVector(dd4hep::rec::Vector3D vec);

TGeoHMatrix ComputeSensorTrafoMatrix(const dd4hep::DDSegmentation::CellID& cellID, const dd4hep::VolumeManager& volumeManager, const TGeoRotation& sensorNormalRotation);

dd4hep::rec::Vector3D GlobalToLocal(const dd4hep::rec::Vector3D& global, const TGeoHMatrix& M);
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