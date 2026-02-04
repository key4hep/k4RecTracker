#pragma once

#include <vector>
#include <numeric>
#include <tuple>

#include "Gaudi/Property.h"

// #include "DDRec/ISurface.h"
// #include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "DDRec/Material.h"

#include "DD4hep/Objects.h"
#include "DD4hep/Volumes.h"
#include "DD4hep/Shapes.h"
#include "DD4hep/DetElement.h"

namespace VTXdigi_tools {

bool ToolTest();

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
dd4hep::rec::Vector3D ComputePixelCenter_Local(const std::array<int, 2> pixelIndex, const dd4hep::rec::ISurface& surface,  const std::array<float, 2> pixelPitch, float depletedRegionDepthCenter);

/** @brief Compute the center position of a given pixel (i_u,i_v) in sensor-local coordinates (u,v,0) */
dd4hep::rec::Vector3D ComputePixelCenter_Local(const std::array<int, 2> pixelIndex, const dd4hep::rec::ISurface& surface, const std::array<float, 2> pixelPitch);



/* -- Pixel Charge Matrix -- */

/** @brief Holds the charge collected in a 2d matrix of pixels surrounding a simHit.
 * 
 * @note Starts from default matrix size of 11x11 pixels, expands dynamically. Expanding is computationally expensive.
 */
class PixelChargeMatrix {
  /* Stores the charge deposited in a (size_u x size_v) pixel matrix around a given origin pixel.
    * In case charge is added outside the matrix bounds, the matrix range is expanded in that direction.
    * The size of the matrix is defined via the Gaudi property MaximumClusterSize. */

  const int m_initialSize = 11; // initial size of the matrix in u and v direction
  const int m_overExpansionStep = 3; // number of additional pixels to expand the matrix by, when a charge is added outside the current bounds
  std::vector<float> m_pixelCharge;
  int m_range_u[2], m_range_v[2]; // Inclusive matrix range. -> size = range[1] - range[0] + 1
  /* range CAN extend into negative values or outside of sensor area, this ensures graceful handling of hits outside inditial bounds. These might be discarded later, if outside of sensor. */
  int m_origin[2]; // origin pixel indices

  public:

    PixelChargeMatrix(int i_origin_u, int i_origin_v);

    void Reset();
    
    inline int GetOriginU() const { return m_origin[0]; }
    inline int GetOriginV() const { return m_origin[1]; }
    inline int GetRangeMin_u() const { return m_range_u[0]; }
    inline int GetRangeMax_u() const { return m_range_u[1]; }
    inline int GetRangeMin_v() const { return m_range_v[0]; }
    inline int GetRangeMax_v() const { return m_range_v[1]; }
    inline int GetSize_u() const { return m_range_u[1] - m_range_u[0] + 1; }
    inline int GetSize_v() const { return m_range_v[1] - m_range_v[0] + 1; }
    inline std::tuple<int, int> GetSize() const;
    inline std::tuple<int, int> GetOrigin() const;

    float GetTotalCharge() const;
    float GetCharge(int i_u, int i_v) const;
    void FillCharge(int i_u, int i_v, float charge);

  private:
    inline int _FindIndex(int i_u, int i_v) const;
    inline bool _OutOfBounds(int i_u, int i_v) const;
    void _ExpandMatrix(int i_u, int i_v);
}; // class PixelChargeMatrix

class SensorChargeMatrix {
  /* Holds the charge collected in a sensor, stored as a map of pixel indices to charge. This is used as an intermediate step before creating the output TrackerHit, which requires a fixed-size matrix. */

  std::vector<int> m_sensorCharge; // map of pixel indices (i_u, i_v) to charge
  std::array<size_t, 2> m_pixelCount; // number of pixels in u and v direction

  public:
    SensorChargeMatrix(std::array<size_t, 2> pixelCount);

    void FillCharge(std::array<int, 2> pixel, int charge);
    int GetCharge(std::array<int, 2> pixel) const;
    int GetTotalCharge() const;
    void Reset();

  private:
    inline bool _OutOfBounds(std::array<int, 2> pixel) const;
    inline int _FindIndex(std::array<int, 2> pixel) const;
}; // class SensorChargeMatrix

} // namespace VTXdigi_tools