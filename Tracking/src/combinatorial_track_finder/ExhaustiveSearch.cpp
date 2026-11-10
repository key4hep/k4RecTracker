// See ExhaustiveSearch.h.
#include "ExhaustiveSearch.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace k4rt {

int regionOfCid(int cid, const std::unordered_map<int, int>* cidRegion) {
  if (cidRegion) {
    auto it = cidRegion->find(cid);
    if (it != cidRegion->end())
      return it->second;
  }
  return regionOfKey(cid);
}

bool RoadCuts::passesRoadCuts(const std::vector<size_t>& gIdxVec) const {
  const int nH = static_cast<int>(gIdxVec.size());

  // azimuth phi of every hit (positions in mm; the angle does not depend on the unit)
  std::vector<double> phi(nH);
  for (int i = 0; i < nH; ++i) {
    const auto p = m_hits[gIdxVec[i]].getPosition();
    phi[i] = std::atan2(p[1], p[0]);
  }

  // Near the beam all phi cuts are skipped: low-pT tracks there turn by a large angle in phi
  // between endcap layers. The polar angle of the combination is that of its mean point:
  // theta = atan2(mean r, mean |z|).
  {
    double sumR = 0.0, sumAbsZ = 0.0;
    for (int i = 0; i < nH; ++i) {
      const auto p = m_hits[gIdxVec[i]].getPosition();
      sumR += std::sqrt(p[0] * p[0] + p[1] * p[1]);
      sumAbsZ += std::abs(p[2]);
    }
    const double thetaCombo = std::atan2(sumR / nH, sumAbsZ / nH);
    const double nearBeamLo = m_config.nearBeamThetaMinDeg * M_PI / 180.0;
    const double nearBeamHi = m_config.nearBeamThetaMaxDeg * M_PI / 180.0;
    if (thetaCombo < nearBeamLo || thetaCombo > nearBeamHi)
      return true;
  }

  // Cut 1: |dphi| between hits of consecutive layers. A pair across the barrel <-> endcap
  // transition is skipped (a large dphi is expected there).
  const double maxConsec = m_config.maxConsecDeltaPhi;
  if (maxConsec > 0.0) {
    for (int i = 0; i < nH - 1; ++i) {
      // layers of the two hits (compositeID 0 if a hit has no entry)
      int cidA = 0, cidB = 0;
      {
        auto itA = m_compositeOfHit.find(gIdxVec[i]);
        auto itB = m_compositeOfHit.find(gIdxVec[i + 1]);
        if (itA != m_compositeOfHit.end())
          cidA = itA->second;
        if (itB != m_compositeOfHit.end())
          cidB = itB->second;
      }
      const bool aBarrel = (regionOfCid(cidA, m_cidRegion) == 0);
      const bool bBarrel = (regionOfCid(cidB, m_cidRegion) == 0);
      if (aBarrel != bBarrel)
        continue;

      // the azimuthal step, wrapped to (-pi, pi] so that steps across phi = +-pi are small
      double dPhi = phi[i + 1] - phi[i];
      dPhi = wrapAngle(dPhi);
      dPhi = std::abs(dPhi);
      if (dPhi > maxConsec)
        return false;
    }
  }

  // Cut 2: total phi spread, max - min of every hit's phi relative to the innermost hit (each
  // difference wrapped to (-pi, pi], so the window may straddle phi = +-pi)
  const double maxSpread = m_config.maxComboPhiSpread;
  if (maxSpread > 0.0) {
    double lo = 0.0, hi = 0.0;
    for (int i = 1; i < nH; ++i) {
      double d = phi[i] - phi[0];
      d = wrapAngle(d);
      if (d < lo)
        lo = d;
      if (d > hi)
        hi = d;
    }
    const double spread = hi - lo;
    if (spread > maxSpread)
      return false;
  }

  // Cut 3: in each endcap, the transverse radius r of a track coming from the IP can only grow
  // (or stay the same) when |z| grows. A hit where r decreases while |z| increases cannot be on
  // the same track (a secondary, or a hit from another track on another disk).
  {
    std::vector<std::pair<double, double>> posEC, negEC; // (|z| mm, r mm) of the hits of each endcap
    for (int i = 0; i < nH; ++i) {
      const auto p = m_hits[gIdxVec[i]].getPosition();
      auto it = m_compositeOfHit.find(gIdxVec[i]);
      if (it == m_compositeOfHit.end())
        continue;
      const int cid = it->second;
      const double r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
      const double absZ = std::abs(p[2]);
      const int region = regionOfCid(cid, m_cidRegion);
      if (region == 1)
        posEC.emplace_back(absZ, r);
      if (region == -1)
        negEC.emplace_back(absZ, r);
    }
    if (!radiusGrowsWithZ(posEC) || !radiusGrowsWithZ(negEC))
      return false;
  }

  return true;
}

bool RoadCuts::radiusGrowsWithZ(std::vector<std::pair<double, double>>& ec) {
  if (ec.size() < 2)
    return true;
  // sorted by |z|, r must never decrease
  std::sort(ec.begin(), ec.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  for (int i = 0; i < (int)ec.size() - 1; ++i) {
    const double dr = ec[i + 1].second - ec[i].second;
    if (dr < 0.0)
      return false;
  }
  return true;
}

double RoadCuts::distanceCm(size_t first, size_t second) const {
  const auto pa = m_hits[first].getPosition();
  const auto pb = m_hits[second].getPosition();
  const double dx = pb[0] - pa[0], dy = pb[1] - pa[1], dz = pb[2] - pa[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz) * 0.1;
}

bool RoadCuts::passesDistanceCuts(const std::vector<size_t>& gIdxVec) const {
  const int nH = static_cast<int>(gIdxVec.size());
  const double maxConsec = m_config.maxConsecHitDist;
  const double maxPair = m_config.maxPairHitDist;

  // consecutive hits (inner -> outer)
  if (maxConsec > 0.0) {
    for (int i = 0; i < nH - 1; ++i) {
      if (distanceCm(gIdxVec[i], gIdxVec[i + 1]) > maxConsec)
        return false;
    }
  }
  // any two hits
  if (maxPair > 0.0) {
    for (int i = 0; i < nH - 1; ++i) {
      for (int j = i + 1; j < nH; ++j) {
        if (distanceCm(gIdxVec[i], gIdxVec[j]) > maxPair)
          return false;
      }
    }
  }
  return true;
}

bool RoadCuts::passesTanLambdaCut(const std::vector<size_t>& gIdxVec) const {
  const double maxDev = m_config.maxComboTanLambdaDev;
  if (maxDev <= 0.0)
    return true;
  const int nH = static_cast<int>(gIdxVec.size());
  if (nH < 3)
    return true;

  // Mixed barrel + endcap combinations have, by geometry, very different dz/dr per pair: the cut
  // is only meaningful for pure barrel or pure endcap segments.
  bool hasBarrel = false, hasEndcap = false;
  for (size_t gi : gIdxVec) {
    auto it = m_compositeOfHit.find(gi);
    if (it == m_compositeOfHit.end())
      continue;
    const int cid = it->second;
    if (regionOfCid(cid, m_cidRegion) == 0)
      hasBarrel = true;
    else
      hasEndcap = true;
    if (hasBarrel && hasEndcap)
      return true;
  }

  // dip of every pair of consecutive hits: tan(lambda) = dz / (transverse distance). For a
  // helix the transverse chord is close to the arc length when the hits are not far apart.
  std::vector<double> tanLPairs;
  tanLPairs.reserve(nH - 1);
  for (int i = 0; i < nH - 1; ++i) {
    const auto pa = m_hits[gIdxVec[i]].getPosition();
    const auto pb = m_hits[gIdxVec[i + 1]].getPosition();
    const double dxy = std::sqrt((pb[0] - pa[0]) * (pb[0] - pa[0]) + (pb[1] - pa[1]) * (pb[1] - pa[1]));
    if (dxy < m_config.minPairTransverseDistanceMm)
      continue; // too close in the transverse plane for a meaningful slope
    tanLPairs.push_back((pb[2] - pa[2]) / dxy);
  }
  if (static_cast<int>(tanLPairs.size()) < 2)
    return true;

  // every pair must be close to the median dip (the median is not pulled by one wrong hit)
  std::vector<double> sorted = tanLPairs;
  std::sort(sorted.begin(), sorted.end());
  const double median = sorted[sorted.size() / 2];
  for (double tl : tanLPairs) {
    if (std::abs(tl - median) > maxDev)
      return false;
  }
  return true;
}

} // namespace k4rt
