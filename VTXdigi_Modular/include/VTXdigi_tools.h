// VTXdigi_Modular/include/VTXdigi_tools.h
#pragma once

#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/RndmGenerators.h"

#include "DDRec/Surface.h"
#include "DDRec/CellIDPositionConverter.h"

#include "edm4hep/SimTrackerHit.h"

#include <queue>
#include <unordered_set>
#include <string_view>
#include <limits>

namespace VTXdigi_tools {

constexpr float kChargePerkeV = 273.97f; // in electrons, for silicon (1 eh-pair ~ 3.65 eV)

/* -- SimHitWrapper -- */

enum class MCParticleLevel {
  Primary, // created in generator
  Secondary, // created in simulation (scattering or decay) but not inside this sensors volume
  Delta // created in simulation, but either (a) MCParticle dropped in ddsim because of low energy or (b) created inside the sensor volume
};

/** @brief Class to contain all information about a simTrackerHit that is needed for the digitization.
 * @note this is where the simTrackerHit is actually stored, everything else (pixelHit / cluster) will store pointers to this. */
class SimHitWrapper {
  edm4hep::SimTrackerHit m_simTrackerHit;
  dd4hep::DDSegmentation::VolumeID m_volumeID; // this is the CellID without segmentation bits
  float m_charge;
  int m_layerNumber;
  mutable dd4hep::rec::Vector3D m_truthPos; // simHit truth position, local coordinates. Mutable because it might be adjusted in const ChargeCollector::FillHit() to account for charge collection biases (see ChargeCollector_impl.h ChargeCollector_LUT::MoveTruthPosition() for more info).
  MCParticleLevel m_mcParticleLevel;

public:
  SimHitWrapper(
    edm4hep::SimTrackerHit simTrackerHit, dd4hep::DDSegmentation::VolumeID volumeID,
    const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder,
    const dd4hep::VolumeManager& volumeManager,
    const std::unique_ptr<dd4hep::rec::CellIDPositionConverter>& cellIDPositionConverter);
  SimHitWrapper(const SimHitWrapper& other) = default;
  SimHitWrapper(SimHitWrapper&& other) = default;
  SimHitWrapper() = default;

  /** @brief Set the truth position of the simHit in local coordinates */
  inline void SetTruthPos(const dd4hep::rec::Vector3D& pos) const { m_truthPos = pos; } // only used for histogramming after filling the hits, so not really a problem that this is mutable

  friend void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept;
  inline const edm4hep::SimTrackerHit* hitPtr() const { return &m_simTrackerHit; }

  /** @brief Access the truth position of the simHit in local coordinates
   * @note might have been adjusted by ChargeCollector::FillHit() to correct for charge collection effects (eg. not completely depleted sensors where charges are only collected close to the upper sensor surface as in TPSCo 65nm CIS). */
  inline const dd4hep::rec::Vector3D truthPos() const { return m_truthPos; }

  inline dd4hep::DDSegmentation::VolumeID volumeID() const { return m_volumeID; }
  inline float charge() const { return m_charge; }
  inline int layer() const { return m_layerNumber; }

  /** @brief Determine whether the simHit was produced by a (a) primary, (b) secondary or (c) delta particle
   * primary - generator level particle
   * secondary - particle produced in Geant4 simulation but not inside this sensors volume
   * delta - particle produced in Geant4 simulation and either (a) dropped in ddsim or (b) produced inside the sensor volume */
  inline MCParticleLevel mcParticleLevel() const { return m_mcParticleLevel; }
};

void swap(SimHitWrapper& a, SimHitWrapper& b) noexcept;

/** @brief Find simHits from different MCParticles. If two simHits originate from the same MCParticle (eg. have shared parents), return only the one from the MCParticle further up the family tree. */
// std::unordered_set<const VTXdigi_tools::SimHitWrapper*>  FindSimHitsWithIndividualParents(const std::vector<SimHitWrapper>& simHits);

/* -- Pixel -- */

/** @brief A pixel in the hit map. Can have multiple contributing simHits */
struct Pixel {
  float charge;
  std::unordered_set<const SimHitWrapper*> simHits;
  std::array<int, 2> index; // This info is saved in (a) the map key, and (b) here inside the Pixel object. This is inefficient. But it makes the code a bit nicer not having to pass the index around separately.

  Pixel(std::array<int, 2> pix) : charge(0.f), index(pix) {
    simHits.reserve(2); // avoid too many reallocations, will rarely see more than 2 simTrackerHits contributing to the same pixel
  }
  Pixel() : charge(0.f), index({-1, -1}) {
    simHits.reserve(2); // avoid too many reallocations, will rarely see more than 2 simTrackerHits contributing to the same pixel
  }
};

/* -- Cluster -- */

/** @brief Holds clusters: charge, pointers to all contributing pixels, pointers to all contributing simHits */
struct Cluster {
  std::vector<const Pixel*> pixels; // raw pointer into HitMap::Pixels. Valid as long as the HitMap is not modified after clustering. See HitMap::ComputeClusters() for details.
  std::unordered_set<const SimHitWrapper*> simHits;
  float charge = 0.f;

  inline int GetSize() const { return pixels.size(); };
  int GetSize(const int axis) const; // axis = 0 for u, 1 for v.

  /** @brief get the charge in the seed pixel (eg. pixel with the highest charge) */
  float GetSeedPixelCharge() const;

  /** @brief Compute the charge-weighted centre-of-gravity of a cluster in terms of pixel inx coordinates
   * @returns CoG in terms of pixel index coordinates */
  std::array<float, 2> ComputeCoG(const bool clusterizeEndPixelsOnly) const;

  /** @brief Compute the uncertainty of the charge-weighted centre-of-gravity of a cluster in terms of pixel index coordinates
   * @param pos This cluster's CoG (to avoid re-computing it) in terms of pixel index coordinates */
  std::array<float, 2> ComputeCoGUncertainty(const std::array<float, 2>& pos) const;
};

/** @brief Get the indices of all direct neighbors of a pixel */
std::array<std::array<int, 2>, 4> GetDirectNeighbors(const std::array<int, 2>& i_uv);
std::array<std::array<int, 2>, 8> GetNeighbors(const std::array<int, 2>& i_uv);

/* -- HitMap -- */

struct Hash_PixIndex {
  size_t operator()(const std::array<int, 2>& i_uv) const noexcept {
    return (static_cast<uint64_t>(i_uv[0]) << 32) ^ static_cast<uint32_t>(i_uv[1]);
  }
};

using PixelMap = std::unordered_map<std::array<int, 2>, Pixel, Hash_PixIndex>;

/** @brief HitMap of all pixel hits on a sensor
 * @note uses a std::unordered_map to only store pixels that have charge, which is more memory efficient for large pixel counts and low occupancy. */
class HitMap {
  /* I tried implementing a vector that contains every pixel, but for large pixel counts this is very memory-inefficient. Instead, this class uses a std::map to only store pixels that have charge. This is more memory efficient for sparse hits, with O(1) simHit/sensor/event this is at least a factor 100 faster that the vector approach. Maybe not the case to ttbar run with O(0.1%) pixel occupancy. */

  PixelMap m_pixels; // hit pixels, stored by value
  std::array<size_t, 2> m_pixCount; // size of the sensor in pixels

public:
  HitMap(std::array<size_t, 2> pixCount);

  /** @brief Add charge and a simHit to a pixel */
  void FillCharge(std::array<int, 2> i_uv, float charge, const SimHitWrapper& simHitWrapper);

  /** @brief For each pixel with charge, vary the charge by an amount drawn from the supplied random generator */
  void ApplyChargeSmearing(const Rndm::Numbers& rndm_charge);

  /** @brief Erase pixels below threshold. If rndm_threshold is given, a per-pixel Gaussian
   *  dispersion is drawn from it and added to the threshold before comparison. */
  void ApplyThreshold(const float threshold, const Rndm::Numbers* rndm_threshold = nullptr);

  /** @brief Get one pixel's collected charge */
  float GetCharge(std::array<int, 2> i_uv) const;

  /** @brief Get the total charge across all pixels */
  float GetTotalCharge() const;

  /** @brief Return a const reference to the internal pixel map */
  inline const PixelMap& Hits() const { return m_pixels; };

  /** @brief Return the number of pixels with charge */
  inline int GetTotalPixelsWithCharge() const { return m_pixels.size(); };

  /** @brief Clusterize hit pixels in the HitMap, using direct neighbors
   * @note Returns Cluster objects whose `pixels` members hold raw pointers into m_pixels. This is safe because ComputeClusters() is const and called only after all FillCharge() insertions are complete -> the map will not rehash while the returned Clusters are alive. Do not call FillCharge() on this HitMap after calling ComputeClusters(). */
  std::vector<Cluster> ComputeClusters() const;

  /** @brief Return a vector of clusters, where each cluster is a single pixel */
  std::vector<Cluster> ComputeClusters_singePixels() const;

  inline void Reset() { m_pixels.clear(); };

private:
  /** @brief Returns true if the pixel is out of bounds */
  inline bool _OutOfBounds(std::array<int, 2> i_uv) const;
}; // class HitMap

/* -- helpers -- */

/** @brief Convert a edm4hep::Vector3d to dd4hep::rec::Vector3D */
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3d vec);
/** @brief Convert a edm4hep::Vector3f to dd4hep::rec::Vector3D */
dd4hep::rec::Vector3D ConvertVector(edm4hep::Vector3f vec);
/** @brief Convert a dd4hep::rec::Vector3D to edm4hep::Vector3d */
edm4hep::Vector3d ConvertVector(dd4hep::rec::Vector3D vec);


/** @brief Compute the transformation matrix from global detector to local sensor frame for a given sensor volume (defined by its volumeID) */
TGeoHMatrix ComputeSensorTrafoMatrix(const dd4hep::DDSegmentation::VolumeID& volumeID, const dd4hep::VolumeManager& volumeManager, const TGeoRotation& sensorNormalRotation);

/** @brief Transform a position from global detector coordinates to sensor-local coordinates, using the sensor transformation matrix */
dd4hep::rec::Vector3D Trafo_global_local(const dd4hep::rec::Vector3D& global, const TGeoHMatrix& M);
/** @brief Transform a position from sensor-local coordinates to global detector coordinates, using the sensor transformation matrix */
dd4hep::rec::Vector3D Trafo_local_global(const dd4hep::rec::Vector3D& local, const TGeoHMatrix& M);

/** @brief Transform a position from sensor-local sensor coordinates to pixel index coordinates
 * @note the origin lies at the centre of the (0,0) pixel, pixel centres are at multiples of one */
std::array<float, 2> Trafo_local_pixIndexCoords(const dd4hep::rec::Vector3D& local, const std::array<float, 2> pixelPitch, const std::array<size_t, 2> pixelCount);

/** @brief Transform a position from pixel index coordinates to sensor-local coordinates */
dd4hep::rec::Vector3D Trafo_pixIndexCoords_local(const std::array<float, 2>& pixIndexCoords, const float w, const std::array<float, 2> pixelPitch, const std::array<size_t, 2> pixelCount);

/** @brief Compute the indices of the pixel (i_u, i_v) that a given sensor-local position lies in */
std::array<int, 2> Trafo_local_pixIndex(const dd4hep::rec::Vector3D& pos, const std::array<float, 2> pixelPitch, const std::array<size_t, 2> pixelCount);

/** @brief Compute the inices of the in-pixel bin (j_u, j_v, j_w) that a given sensor-local position lies in */
std::array<int, 3> Trafo_local_inpixIndex(const dd4hep::rec::Vector3D& pos, const std::array<int, 3>& binCount, const std::array<float, 2>& pixelPitch, const std::array<float, 3>& activeVolumeDimensions);

/** @brief Transform a position from pixel index coordinates to sensor-local coordinates
 * @note The w coordinate is set to depletedRegionDepthCenter.*/
dd4hep::rec::Vector3D Trafo_pixIndex_local(const std::array<int, 2> pixelIndex, const std::array<float, 2> sensorLength,  const std::array<float, 2> pixelPitch, float depletedRegionDepthCenter);
/** @brief Transform a position from pixel index coordinates to sensor-local coordinates, setting w=0 */
dd4hep::rec::Vector3D Trafo_pixIndex_local(const std::array<int, 2> pixelIndex, const std::array<float, 2> sensorLength, const std::array<float, 2> pixelPitch);

/** @brief Transform a position from pixel index coordinates to sensor-local coordinates
 * @note The w coordinate is set to depletedRegionDepthCenter */
dd4hep::rec::Vector3D Trafo_pixIndex_local(const std::array<float, 2> index, const std::array<float, 2> sensorLength,  const std::array<float, 2> pixelPitch, float depletedRegionDepthCenter);
/** @brief Transform a position from pixel index coordinates to sensor-local coordinates, setting w=0 */
dd4hep::rec::Vector3D Trafo_pixIndex_local(const std::array<float, 2> index, const std::array<float, 2> sensorLength, const std::array<float, 2> pixelPitch);


int GetLayer(const dd4hep::DDSegmentation::VolumeID& volumeID, const std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder>& cellIdDecoder);

/* -- Binning tools -- */

/** @brief Given a histogram definition (x0, binWidth, nBins) and a value x, compute the bin index i in which x falls.
 * @return Int, -1 if x is out of range.
 * @note Bins are 0-indexed (vs ROOT's 1-indexing) */
int ComputeBinIndex(float x, float binX0, float binWidth, int binN);

/** @brief Given a binning definition (x0, binWidth, nBins), compute the center position of a given bin index i in the histogram
 * @return Float, the center position of the bin
 * @note Bins are 0-indexed (vs ROOT's 1-indexing) */
float ComputeBinCenter(int i, float binX0, float binWidth);
/** @brief Given a binning definition (x0, x1, nBins), compute the center position of a given bin index i in the histogram
 * @return Float, the center position of the bin
 * @note Bins are 0-indexed (vs ROOT's 1-indexing) */
float ComputeBinCenter(int i, float binX0, float binX1, int binN);
} // namespace VTXdigi_tools
