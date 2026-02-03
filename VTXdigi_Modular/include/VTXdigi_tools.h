#pragma once

#include <vector>
#include <numeric>
#include <tuple>

namespace VTXdigi_tools {

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

} // namespace VTXdigi_tools