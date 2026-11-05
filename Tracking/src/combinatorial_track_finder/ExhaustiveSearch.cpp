// See ExhaustiveSearch.h.
#include "ExhaustiveSearch.h"

#include "GaudiKernel/MsgStream.h"

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

bool passesRoadCuts(const std::vector<size_t>& gIdxVec, const RoadCutContext& ctx) {
  const auto& hits = ctx.hits;
  const auto& gIdxToComposite = ctx.gIdxToComposite;
  const int nH = static_cast<int>(gIdxVec.size());

  // Pre-compute phi for each hit (mm positions, ratio is unit-independent)
  std::vector<double> phi(nH);
  for (int i = 0; i < nH; ++i) {
    const auto p = hits[gIdxVec[i]].getPosition();
    phi[i] = std::atan2(p[1], p[0]);
  }

  // Cut 1: consecutive-layer |dphi|
  // Bypass conditions:
  //   - near beam pipe: combo theta < 15 deg or > 160 deg -> bypass ALL phi road cuts
  //     (low-pT tracks here spiral in phi between endcap layers, producing large dphi)
  //   - barrel<->endcap boundary pair -> skip that pair only (transition region)
  {
    double sumR = 0.0, sumAbsZ = 0.0;
    for (int i = 0; i < nH; ++i) {
      const auto p = hits[gIdxVec[i]].getPosition();
      sumR += std::sqrt(p[0] * p[0] + p[1] * p[1]);
      sumAbsZ += std::abs(p[2]);
    }
    const double thetaCombo = std::atan2(sumR / nH, sumAbsZ / nH);
    constexpr double kNearBeamLo = 15.0 * M_PI / 180.0;
    constexpr double kNearBeamHi = 160.0 * M_PI / 180.0;
    if (thetaCombo < kNearBeamLo || thetaCombo > kNearBeamHi)
      return true;
  }
  const double maxConsec = ctx.cfg.maxConsecDeltaPhi;
  if (maxConsec > 0.0) {
    for (int i = 0; i < nH - 1; ++i) {
      int cidA = 0, cidB = 0;
      {
        auto itA = gIdxToComposite.find(gIdxVec[i]);
        auto itB = gIdxToComposite.find(gIdxVec[i + 1]);
        if (itA != gIdxToComposite.end())
          cidA = itA->second;
        if (itB != gIdxToComposite.end())
          cidB = itB->second;
      }
      const bool aBarrel = (regionOfCid(cidA, ctx.cidRegion) == 0);
      const bool bBarrel = (regionOfCid(cidB, ctx.cidRegion) == 0);
      // Barrel<->endcap transition: large dphi is expected - skip this pair
      if (aBarrel != bBarrel)
        continue;

      double dPhi = phi[i + 1] - phi[i];
      dPhi = wrapAngle(dPhi);
      dPhi = std::abs(dPhi);
      if (dPhi > maxConsec) {
        if (ctx.dbg)
          *ctx.dbg << "    road FAIL (consec dphi): gIdx=" << gIdxVec[i] << "->" << gIdxVec[i + 1] << "  dPhi=" << dPhi
                   << " rad  cut=" << maxConsec << " rad" << endmsg;
        return false;
      }
    }
  }

  // Cut 2: total phi spread (max - min relative to innermost hit, wrapped to [-pi,pi])
  const double maxSpread = ctx.cfg.maxComboPhiSpread;
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
    double spread = hi - lo;
    if (spread > maxSpread) {
      if (ctx.dbg)
        *ctx.dbg << "    road FAIL (phi spread): spread=" << spread << " rad  cut=" << maxSpread << " rad" << endmsg;
      return false;
    }
  }

  // Cut 3: endcap z/r monotonicity.
  // For a forward-going track the transverse radius r must increase (or stay constant)
  // as |z| increases. A hit where r decreases while |z| increases is geometrically
  // inconsistent with any single forward track and signals a secondary or noise hit
  // merged from a different disk.
  {
    std::vector<std::pair<double, double>> posEC, negEC; // (|z| mm, r mm)
    for (int i = 0; i < nH; ++i) {
      const auto p = hits[gIdxVec[i]].getPosition();
      auto it = gIdxToComposite.find(gIdxVec[i]);
      if (it == gIdxToComposite.end())
        continue;
      const int cid = it->second;
      const double r = std::sqrt(p[0] * p[0] + p[1] * p[1]);
      const double absZ = std::abs(p[2]);
      const int region = regionOfCid(cid, ctx.cidRegion);
      if (region == 1)
        posEC.emplace_back(absZ, r);
      if (region == -1)
        negEC.emplace_back(absZ, r);
    }
    auto checkZRMono = [&](std::vector<std::pair<double, double>>& ec) -> bool {
      if (ec.size() < 2)
        return true;
      std::sort(ec.begin(), ec.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
      for (int i = 0; i < (int)ec.size() - 1; ++i) {
        const double dr = ec[i + 1].second - ec[i].second;
        if (dr < 0.0) {
          if (ctx.dbg)
            *ctx.dbg << "    road FAIL (endcap z/r mono): dr=" << dr << " mm  |z|=[" << ec[i].first << "->"
                     << ec[i + 1].first << "] mm" << endmsg;
          return false;
        }
      }
      return true;
    };
    if (!checkZRMono(posEC) || !checkZRMono(negEC))
      return false;
  }

  return true;
}

bool passesDistanceCuts(const std::vector<size_t>& gIdxVec, const RoadCutContext& ctx) {
  const auto& hits = ctx.hits;
  const int nH = static_cast<int>(gIdxVec.size());
  const double maxConsec = ctx.cfg.maxConsecHitDist;
  const double maxPair = ctx.cfg.maxPairHitDist;

  if (maxConsec > 0.0) {
    for (int i = 0; i < nH - 1; ++i) {
      const auto pa = hits[gIdxVec[i]].getPosition();
      const auto pb = hits[gIdxVec[i + 1]].getPosition();
      double dx = pb[0] - pa[0], dy = pb[1] - pa[1], dz = pb[2] - pa[2];
      double d = std::sqrt(dx * dx + dy * dy + dz * dz) * 0.1; // mm -> cm
      if (d > maxConsec) {
        if (ctx.dbg)
          *ctx.dbg << "    dist FAIL (consec): hit " << gIdxVec[i] << "->" << gIdxVec[i + 1] << "  d=" << d
                   << " cm  cut=" << maxConsec << " cm" << endmsg;
        return false;
      }
    }
  }
  if (maxPair > 0.0) {
    for (int i = 0; i < nH - 1; ++i) {
      for (int j = i + 1; j < nH; ++j) {
        const auto pa = hits[gIdxVec[i]].getPosition();
        const auto pb = hits[gIdxVec[j]].getPosition();
        double dx = pb[0] - pa[0], dy = pb[1] - pa[1], dz = pb[2] - pa[2];
        double d = std::sqrt(dx * dx + dy * dy + dz * dz) * 0.1; // mm -> cm
        if (d > maxPair) {
          if (ctx.dbg)
            *ctx.dbg << "    dist FAIL (pair): hit " << gIdxVec[i] << "<->" << gIdxVec[j] << "  d=" << d
                     << " cm  cut=" << maxPair << " cm" << endmsg;
          return false;
        }
      }
    }
  }
  return true;
}

bool passesTanLambdaCut(const std::vector<size_t>& gIdxVec, const RoadCutContext& ctx) {
  const auto& hits = ctx.hits;
  const auto& gIdxToComposite = ctx.gIdxToComposite;
  const double maxDev = ctx.cfg.maxComboTanLambdaDev;
  if (maxDev <= 0.0)
    return true;
  const int nH = static_cast<int>(gIdxVec.size());
  if (nH < 3)
    return true;

  // Mixed barrel+endcap combos have very different tanLambda per pair by geometry;
  // the cut is only meaningful for pure-barrel or pure-endcap track segments.
  bool hasBarrel = false, hasEndcap = false;
  for (size_t gi : gIdxVec) {
    auto it = gIdxToComposite.find(gi);
    if (it == gIdxToComposite.end())
      continue;
    int cid = it->second;
    if (regionOfCid(cid, ctx.cidRegion) == 0)
      hasBarrel = true;
    else
      hasEndcap = true;
    if (hasBarrel && hasEndcap)
      return true; // mixed combo - skip cut
  }

  std::vector<double> tanLPairs;
  tanLPairs.reserve(nH - 1);
  for (int i = 0; i < nH - 1; ++i) {
    const auto pa = hits[gIdxVec[i]].getPosition();
    const auto pb = hits[gIdxVec[i + 1]].getPosition();
    double dxy = std::sqrt((pb[0] - pa[0]) * (pb[0] - pa[0]) + (pb[1] - pa[1]) * (pb[1] - pa[1]));
    if (dxy < 1.0)
      continue; // skip degenerate pairs (< 1 mm transverse separation)
    tanLPairs.push_back((pb[2] - pa[2]) / dxy);
  }
  if (static_cast<int>(tanLPairs.size()) < 2)
    return true;

  std::vector<double> sorted = tanLPairs;
  std::sort(sorted.begin(), sorted.end());
  double median = sorted[sorted.size() / 2];

  for (double tl : tanLPairs) {
    if (std::abs(tl - median) > maxDev) {
      if (ctx.dbg)
        *ctx.dbg << "    tanLambda FAIL: pair tanL=" << tl << " median=" << median << " |dev|=" << std::abs(tl - median)
                 << " cut=" << maxDev << endmsg;
      return false;
    }
  }
  return true;
}

} // namespace k4rt
