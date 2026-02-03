#include "VTXdigi_tools.h"
#include <iostream>

namespace VTXdigi_tools {

PixelChargeMatrix::PixelChargeMatrix(int i_origin_u, int i_origin_v) : m_pixelCharge(m_initialSize * m_initialSize, 0.f), m_origin{ i_origin_u, i_origin_v } {
    // At first, the range is centered around the origin
  const int halfSize = (m_initialSize - 1) / 2;
  m_range_u[0] = m_origin[0] - halfSize;
  m_range_u[1] = m_origin[0] + halfSize;
  m_range_v[0] = m_origin[1] - halfSize;
  m_range_v[1] = m_origin[1] + halfSize;
}

void PixelChargeMatrix::Reset() {
  std::fill(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f);
}

inline std::tuple<int, int> PixelChargeMatrix::GetSize() const {
  return std::make_tuple(GetSize_u(), GetSize_v());
}
inline std::tuple<int, int> PixelChargeMatrix::GetOrigin() const {
  return std::make_tuple(m_origin[0], m_origin[1]);
}

float PixelChargeMatrix::GetTotalCharge() const {
  return std::accumulate(m_pixelCharge.begin(), m_pixelCharge.end(), 0.f);
}

float PixelChargeMatrix::GetCharge(int i_u, int i_v) const {
  if (_OutOfBounds(i_u, i_v)) {
    throw std::runtime_error("PixelChargeMatrix::GetCharge: pixel i_u or i_v ( " + std::to_string(i_u) + ", " + std::to_string(i_v) + ") out of range");
  }
  return m_pixelCharge[_FindIndex(i_u, i_v)];
}

void PixelChargeMatrix::FillCharge(int i_u, int i_v, float charge) {
  if (_OutOfBounds(i_u, i_v)) {
    _ExpandMatrix(i_u, i_v);
  }
  m_pixelCharge[_FindIndex(i_u, i_v)] += charge;
}

inline int PixelChargeMatrix::_FindIndex(int i_u, int i_v) const {
  int i_u_rel = i_u - m_range_u[0];
  int i_v_rel = i_v - m_range_v[0];
  return i_u_rel + i_v_rel * GetSize_u();
}

inline bool PixelChargeMatrix::_OutOfBounds(int i_u, int i_v) const {
  return (
    i_u < m_range_u[0]
    || i_u > m_range_u[1]
    || i_v < m_range_v[0]
    || i_v > m_range_v[1]);
}

void PixelChargeMatrix::_ExpandMatrix(int i_u, int i_v) {
  /* Expand the matrix to include (i_u, i_v) and an excess of m_expansionStep pixels. */

  int rangeNew_u[2] = { m_range_u[0], m_range_u[1] }; 
  int rangeNew_v[2] = { m_range_v[0], m_range_v[1] };

  // TODO: optimise the expansion logic below, to expand to two directions at once if hit is close to corner ~ Jona 2025-10

  /* Expand it in each direction if necessary */
  if (i_u < m_range_u[0]) {
    rangeNew_u[0] = i_u - m_overExpansionStep;
  } else if (i_u > m_range_u[1]) {
    rangeNew_u[1] = i_u + m_overExpansionStep;
  }
  if (i_v < m_range_v[0]) {
    rangeNew_v[0] = i_v - m_overExpansionStep;
  } else if (i_v > m_range_v[1]) {
    rangeNew_v[1] = i_v + m_overExpansionStep;
  }

  const int newSizeU = rangeNew_u[1] - rangeNew_u[0] + 1;
  const int newSizeV = rangeNew_v[1] - rangeNew_v[0] + 1;

  std::vector<float> newPixelCharge(newSizeU * newSizeV, 0.f);

  /* Copy old charges into new vector, row by row */
  for (int row = 0; row < GetSize_v(); ++row) {
    std::copy(
      m_pixelCharge.begin() + row * GetSize_u(),
      m_pixelCharge.begin() + (row + 1) * GetSize_u(),
      newPixelCharge.begin() + (row + (m_range_v[0] - rangeNew_v[0])) * newSizeU + (m_range_u[0] - rangeNew_u[0])
    );
  }

  m_pixelCharge.swap(newPixelCharge);
  std::copy(rangeNew_u, rangeNew_u + 2, m_range_u);
  std::copy(rangeNew_v, rangeNew_v + 2, m_range_v);
} // _ExpandMatrix()

} // namespace VTXdigi_tools


