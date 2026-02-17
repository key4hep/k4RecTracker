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



namespace VTXdigi_tools {



/* -- Tool tests -- */

bool ToolTest();


/* -- SimHit -- */

class Hit {
  /* any member that is added needs to be added to swap, too! */
  edm4hep::SimTrackerHit m_simHit;
  dd4hep::rec::ISurface* m_surface;
  
  dd4hep::DDSegmentation::CellID m_cellID; // cellID (without segmentation bits)
  int m_charge;
  int m_layerNumber;
  int m_nSegments = 0;


public:
  Hit(edm4hep::SimTrackerHit simHit, const dd4hep::rec::SurfaceMap* surfaceMap, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder);
  Hit(const Hit& other) = default;
  Hit(Hit&& other) = default;

  friend void swap(Hit& a, Hit& b) noexcept;

  inline const edm4hep::SimTrackerHit& simHit() const { return m_simHit; }
  inline dd4hep::rec::ISurface* surface() const { return m_surface; }

  inline dd4hep::DDSegmentation::CellID cellID() const { return m_cellID; }
  inline int charge() const { return m_charge; }
  inline int layer() const { return m_layerNumber; }
  inline int nSegments() const { return m_nSegments; }
};

void swap(Hit& a, Hit& b) noexcept;

/* -- helpers -- */

void CreateDigiHit(const edm4hep::SimTrackerHit& simHit, edm4hep::TrackerHitPlaneCollection& digiHits, edm4hep::TrackerHitSimTrackerHitLinkCollection& digiHitLinks, const dd4hep::rec::Vector3D& position, const int charge);

dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec);
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3f vec);
edm4hep::Vector3d ConvertVector(dd4hep::rec::Vector3D vec);

TGeoHMatrix ComputeSensorTrafoMatrix(const dd4hep::DDSegmentation::CellID& cellID, const dd4hep::VolumeManager& volumeManager, const TGeoRotation& sensorNormalRotation);

dd4hep::rec::Vector3D GlobalToLocal(const dd4hep::rec::Vector3D& global, const TGeoHMatrix& M);
dd4hep::rec::Vector3D LocalToGlobal(const dd4hep::rec::Vector3D& local, const TGeoHMatrix& M);

dd4hep::DDSegmentation::CellID GetCellID_short(const edm4hep::SimTrackerHit& simHit);

int GetLayer(const dd4hep::DDSegmentation::CellID& cellID, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder);
int GetLayer(const edm4hep::SimTrackerHit& simHit, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder);


/* -- Binning tools -- */

/** @brief Given a histogram definition (x0, binWidth, nBins) and a value x, compute the bin index i in which x falls.
 * @return Int, -1 if x is out of range.
 * @note Bins are 0-indexed (vs ROOT's 1-indexing) */
int ComputeBinIndex(float x, float binX0, float binWidth, int binN);

/** @brief Compute the pixel indices (i_u, i_v) for a given (local) position inside the sensor */
std::array<int, 2> ComputePixelIndices(const dd4hep::rec::Vector3D& pos, const std::array<float, 2> pixelPitch, const std::array<size_t, 2> pixelCount);

/** @brief Compute the in-pixel indices (j_u, j_v, j_w) for a given (local) position inside the pixel and layer index
 *  @note Assumption: each layer has only 1 type if sensor */
std::array<int, 3> ComputeInPixelIndices(const dd4hep::rec::Vector3D& pos, const std::array<size_t, 3> binCount, const std::array<float, 2> pixelPitch, const std::array<float, 2> sensorLength, const float sensorThickness);

/** @brief Compute the center position of a given pixel (i_u,i_v) in sensor-local coordinates (u,v,w) 
 * 
 * @note The w coordinate is set to depletedRegionDepthCenter. 0 for center, +25 for sensor surface, +20 for TPSCo 65nm maps.
*/
dd4hep::rec::Vector3D ComputePixelPos_local(const std::array<int, 2> pixelIndex, const std::array<float, 2> sensorLength,  const std::array<float, 2> pixelPitch, float depletedRegionDepthCenter);
/** @brief Compute the center position of a given pixel (i_u,i_v) in sensor-local coordinates (u,v,0) */
dd4hep::rec::Vector3D ComputePixelPos_local(const std::array<int, 2> pixelIndex, const std::array<float, 2> sensorLength, const std::array<float, 2> pixelPitch);


/** @brief Compute the position of a given (pixel-)index (i_u,i_v) in sensor-local coordinates (u,v,0) .
 * @note index 0 indicates the center of the pixel, index -0.5 the lower edge and +0.5 the upper edge.  
 * @note Does not check if the position is within the sensor bounds!
 * @note The w coordinate is set to depletedRegionDepthCenter. 0 for center, +25 for sensor surface, +20 for TPSCo 65nm maps.
 */
dd4hep::rec::Vector3D ComputePos_local(const std::array<float, 2> index, const std::array<float, 2> sensorLength,  const std::array<float, 2> pixelPitch, float depletedRegionDepthCenter);
/** @brief Compute the position of a given (pixel-)index (i_u,i_v) in sensor-local coordinates (u,v,0) .
 * @note index 0 indicates the center of the pixel, index -0.5 the lower edge and +0.5 the upper edge.  
 * @note Does not check if the position is within the sensor bounds!
 */
dd4hep::rec::Vector3D ComputePos_local(const std::array<float, 2> index, const std::array<float, 2> sensorLength, const std::array<float, 2> pixelPitch);


struct PixelHit {
  std::array<int, 2> index; // pixel indices i_u, i_v
  int charge;
};

class HitMap {
  /* I tried implementing a vector that contains every pixel, but for large pixel counts this is very memory-inefficient. Instead, this class uses a std::map to only store pixels that have charge. This is more memory efficient for sparse hits, with O(1) simHit/sensor/event this is at least a factor 100 faster that the vector approach. Maybe not the case to ttbar run with O(0.1%) pixel occupancy. */
  std::unordered_map< int, int> m_pixCharge; // map of pixel indices (i_u, i_v) to charge. std::array<int,2> is not hashable, so we need an index
  std::array<size_t, 2> m_pixCount; // number of pixels in u and v direction

  /* TODO: check if int is enough for pixel indices. 65k could be too small -> use size_t*/

public:
  HitMap(std::array<size_t, 2> pixCount);

  void FillCharge(std::array<int, 2> pix, int charge);
  int GetCharge(std::array<int, 2> pix) const;
  int GetTotalCharge() const;
  std::vector<PixelHit> GetPixelsWithCharge() const;
  inline int GetTotalPixelsWithCharge() const;
  inline void Reset();

private:
  inline bool _OutOfBounds(std::array<int, 2> pix) const;
  inline int _Index(std::array<int, 2> pix) const;
  inline std::array<int, 2> _Index(int index) const;
}; // class HitMap

} // namespace VTXdigi_tools