// See ExhaustiveSearch.h.
#include "ExhaustiveSearch.h"
#include "LayerGeometry.h"
#include "TrackingMath.h"

#include <algorithm>
#include <cmath>

namespace k4rt {

int regionOfLayer(int layerKey, const std::unordered_map<int, int>* regionOfLayerKey) {
  if (regionOfLayerKey) {
    auto entry = regionOfLayerKey->find(layerKey);
    if (entry != regionOfLayerKey->end())
      return entry->second;
  }
  return regionFromKey(layerKey);
}

bool RoadCuts::passesRoadCuts(const std::vector<size_t>& hitIndices) const {
  const int numberOfHits = static_cast<int>(hitIndices.size());

  // near the beam (polar angle of the mean point, atan2(mean r, mean |z|)) the phi cuts are skipped
  double sumR = 0.0;
  double sumAbsZ = 0.0;
  for (size_t hitIndex : hitIndices) {
    sumR += m_hitGeometry[hitIndex].transverseRadius;
    sumAbsZ += m_hitGeometry[hitIndex].absoluteZ;
  }
  const double thetaCombination = std::atan2(sumR / numberOfHits, sumAbsZ / numberOfHits);
  const double nearBeamLow = m_config.nearBeamThetaMinDeg * M_PI / 180.0;
  const double nearBeamHigh = m_config.nearBeamThetaMaxDeg * M_PI / 180.0;
  if (thetaCombination < nearBeamLow || thetaCombination > nearBeamHigh)
    return true;

  // cut 1: |dphi| of consecutive layers (wrapped to (-pi, pi]); barrel <-> endcap pairs skipped
  if (m_config.maxConsecDeltaPhi > 0.0) {
    for (int i = 0; i + 1 < numberOfHits; ++i) {
      const HitGeometry& inner = m_hitGeometry[hitIndices[i]];
      const HitGeometry& outer = m_hitGeometry[hitIndices[i + 1]];
      if ((inner.region == 0) != (outer.region == 0))
        continue;
      const double deltaPhi = std::abs(wrapAngle(outer.azimuth - inner.azimuth));
      if (deltaPhi > m_config.maxConsecDeltaPhi)
        return false;
    }
  }

  // cut 2: phi spread, max - min of every phi relative to the innermost hit (window may cross +-pi)
  if (m_config.maxComboPhiSpread > 0.0) {
    const double innermostAzimuth = m_hitGeometry[hitIndices[0]].azimuth;
    double lowest = 0.0;
    double highest = 0.0;
    for (int i = 1; i < numberOfHits; ++i) {
      const double difference = wrapAngle(m_hitGeometry[hitIndices[i]].azimuth - innermostAzimuth);
      lowest = std::min(lowest, difference);
      highest = std::max(highest, difference);
    }
    if (highest - lowest > m_config.maxComboPhiSpread)
      return false;
  }

  // cut 3: in each endcap r can only grow with |z| for a track from the IP
  std::vector<std::pair<double, double>> positiveEndcapHits; // (|z|, r) in mm
  std::vector<std::pair<double, double>> negativeEndcapHits;
  for (size_t hitIndex : hitIndices) {
    const HitGeometry& hit = m_hitGeometry[hitIndex];
    if (hit.region == 1)
      positiveEndcapHits.emplace_back(hit.absoluteZ, hit.transverseRadius);
    if (hit.region == -1)
      negativeEndcapHits.emplace_back(hit.absoluteZ, hit.transverseRadius);
  }
  return radiusGrowsWithZ(positiveEndcapHits) && radiusGrowsWithZ(negativeEndcapHits);
}

bool RoadCuts::radiusGrowsWithZ(std::vector<std::pair<double, double>>& endcapHits) {
  if (endcapHits.size() < 2)
    return true;
  std::sort(endcapHits.begin(), endcapHits.end(),
            [](const auto& left, const auto& right) { return left.first < right.first; });
  for (size_t i = 0; i + 1 < endcapHits.size(); ++i) {
    if (endcapHits[i + 1].second - endcapHits[i].second < 0.0)
      return false;
  }
  return true;
}

double RoadCuts::distanceCm(size_t first, size_t second) const {
  const HitGeometry& firstHit = m_hitGeometry[first];
  const HitGeometry& secondHit = m_hitGeometry[second];
  const double dx = secondHit.x - firstHit.x, dy = secondHit.y - firstHit.y, dz = secondHit.z - firstHit.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz) * kCmPerMm;
}

bool RoadCuts::passesDistanceCuts(const std::vector<size_t>& hitIndices) const {
  const size_t numberOfHits = hitIndices.size();
  if (m_config.maxConsecHitDist > 0.0) { // consecutive hits
    for (size_t i = 0; i + 1 < numberOfHits; ++i) {
      if (distanceCm(hitIndices[i], hitIndices[i + 1]) > m_config.maxConsecHitDist)
        return false;
    }
  }
  if (m_config.maxPairHitDist > 0.0) { // any two hits
    for (size_t i = 0; i + 1 < numberOfHits; ++i) {
      for (size_t j = i + 1; j < numberOfHits; ++j) {
        if (distanceCm(hitIndices[i], hitIndices[j]) > m_config.maxPairHitDist)
          return false;
      }
    }
  }
  return true;
}

bool RoadCuts::passesTanLambdaCut(const std::vector<size_t>& hitIndices) const {
  const double maxDeviation = m_config.maxComboTanLambdaDev;
  const size_t numberOfHits = hitIndices.size();
  if (maxDeviation <= 0.0 || numberOfHits < 3)
    return true;

  // mixed barrel + endcap combinations have very different dz/dr per pair by geometry: they pass
  bool hasBarrel = false;
  bool hasEndcap = false;
  for (size_t hitIndex : hitIndices) {
    if (m_hitGeometry[hitIndex].region == 0)
      hasBarrel = true;
    else
      hasEndcap = true;
  }
  if (hasBarrel && hasEndcap)
    return true;

  // dip of each consecutive pair: dz / transverse distance (close to the arc length for near hits)
  std::vector<double> pairTanLambdas;
  pairTanLambdas.reserve(numberOfHits - 1);
  for (size_t i = 0; i + 1 < numberOfHits; ++i) {
    const HitGeometry& inner = m_hitGeometry[hitIndices[i]];
    const HitGeometry& outer = m_hitGeometry[hitIndices[i + 1]];
    const double transverseDistance =
        std::sqrt((outer.x - inner.x) * (outer.x - inner.x) + (outer.y - inner.y) * (outer.y - inner.y));
    if (transverseDistance < m_config.minPairTransverseDistanceMm)
      continue; // too close in xy for a slope
    pairTanLambdas.push_back((outer.z - inner.z) / transverseDistance);
  }
  if (pairTanLambdas.size() < 2)
    return true;

  // every pair close to the median (not pulled by one wrong hit)
  std::vector<double> sorted = pairTanLambdas;
  std::sort(sorted.begin(), sorted.end());
  const double median = sorted[sorted.size() / 2];
  for (double pairTanLambda : pairTanLambdas) {
    if (std::abs(pairTanLambda - median) > maxDeviation)
      return false;
  }
  return true;
}

} // namespace k4rt
