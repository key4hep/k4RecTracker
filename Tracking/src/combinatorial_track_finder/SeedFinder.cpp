/**
 * SeedFinder: the SeedAndExtend search of CombinatorialTrackFinder.
 * See SeedFinder.h for an overview of the algorithm.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

#include "SeedFinder.h"
#include "TrackingMath.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <unordered_set>

namespace k4rt {

namespace {
  Eigen::Vector3d posCm(const edm4hep::TrackerHitPlane& h) {
    const auto& p = h.getPosition();
    return {p[0] * 0.1, p[1] * 0.1, p[2] * 0.1};
  }
  // All size-n subsets of v[first..first+len) that contain v[first], order preserved.
  void windowSubsets(const std::vector<int>& v, size_t first, size_t len, int n, std::set<std::vector<int>>& out) {
    if (first + n > v.size())
      return;
    const size_t end = std::min(v.size(), first + len);
    std::vector<int> cur{v[first]};
    std::function<void(size_t)> rec = [&](size_t i) {
      if (static_cast<int>(cur.size()) == n) {
        out.insert(cur);
        return;
      }
      for (size_t j = i; j < end; ++j) {
        cur.push_back(v[j]);
        rec(j + 1);
        cur.pop_back();
      }
    };
    rec(first + 1);
  }
  // Closed-form version of the post-fit checks for a 3-hit seed, applied BEFORE any fit.
  // For 3 points the circle fit passes through all of them, so the analytic 3-point circle is
  // the fitted circle and the radius, travel-order and z-residual checks can be evaluated
  // exactly here. A 1% tolerance absorbs rounding; the full checks still run afterwards and
  // remain authoritative, so this only removes work and never changes which seeds survive.
  // Collinear points (no finite circle) are passed through to the full path unchanged.
  bool tripletPrecheck(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, double maxZRes) {
    const double d = 2.0 * (a.x() * (b.y() - c.y()) + b.x() * (c.y() - a.y()) + c.x() * (a.y() - b.y()));
    if (std::abs(d) < 1e-12)
      return true;
    const double a2 = a.squaredNorm() - a.z() * a.z(), b2 = b.squaredNorm() - b.z() * b.z(),
                 c2 = c.squaredNorm() - c.z() * c.z();
    const double ux = (a2 * (b.y() - c.y()) + b2 * (c.y() - a.y()) + c2 * (a.y() - b.y())) / d;
    const double uy = (a2 * (c.x() - b.x()) + b2 * (a.x() - c.x()) + c2 * (b.x() - a.x())) / d;
    const double R = std::hypot(a.x() - ux, a.y() - uy);
    Helix h;
    std::vector<double> s;
    if (!helixFromHits(ux, uy, R, {a, b, c}, h, &s))
      return true;
    if (!(s[1] > s[0]) || !(s[2] > s[1]))
      return false;
    double zmax = 0.0;
    for (int k = 0; k < 3; ++k) {
      const Eigen::Vector3d& p = (k == 0) ? a : (k == 1) ? b : c;
      zmax = std::max(zmax, std::abs(p.z() - (h.zRef + h.tanL * s[k])));
    }
    return zmax <= maxZRes * 1.01 + 1e-6;
  }
  // Does `pn` lie on the helix through a, b, c (ordered along the track), within
  // window + slope * (arc length from c to pn)? Collinear anchors cannot define a circle; the
  // hit is then accepted and left to the seed fit.
  bool hitOnTrajectory(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
                       const Eigen::Vector3d& pn, double window, double slope) {
    const double d = 2.0 * (a.x() * (b.y() - c.y()) + b.x() * (c.y() - a.y()) + c.x() * (a.y() - b.y()));
    if (std::abs(d) < 1e-12)
      return true;
    const double a2 = a.x() * a.x() + a.y() * a.y(), b2 = b.x() * b.x() + b.y() * b.y(),
                 c2 = c.x() * c.x() + c.y() * c.y();
    const double ux = (a2 * (b.y() - c.y()) + b2 * (c.y() - a.y()) + c2 * (a.y() - b.y())) / d;
    const double uy = (a2 * (c.x() - b.x()) + b2 * (a.x() - c.x()) + c2 * (b.x() - a.x())) / d;
    const double R = std::hypot(a.x() - ux, a.y() - uy);
    Helix h;
    std::vector<double> s;
    if (!helixFromHits(ux, uy, R, {a, b, c}, h, &s))
      return true;
    const double sn = arcLengthNear(h, pn, s.back());
    if (!(sn > s.back()))
      return false; // behind the chain: wrong travel order
    const double tol = window + slope * (sn - s.back());
    return (pn - h.at(sn)).norm() <= tol;
  }
} // namespace

SeedFinder::SeedFinder(const LayerGeometry& geo, SeedConfig cfg, FitFn fit, PlaneFn plane)
    : m_geo(geo), m_cfg(std::move(cfg)), m_fit(std::move(fit)), m_plane(std::move(plane)) {}

std::vector<std::string> SeedFinder::prepare() {
  std::vector<std::string> problems;
  m_sets.clear();
  m_setsByN.clear();
  m_passes.clear();
  for (int p : m_cfg.passes) {
    if (p < 3)
      problems.push_back("a seed size < 3 was raised to 3");
    m_passes.push_back(std::max(3, p));
  }
  for (int n : m_passes) {
    if (m_setsByN.count(n))
      continue;
    m_setsByN[n] = buildSets(n, problems);
  }
  {
    // inside-out: order windows by the closest approach to the IP of their layers, in turn
    auto rank = [&](int cid) {
      const auto* d = m_geo.find(cid);
      return d ? d->dMin_cm : 1e9;
    };
    for (auto& [n, sets] : m_setsByN)
      std::stable_sort(sets.begin(), sets.end(), [&](const std::vector<int>& a, const std::vector<int>& b) {
        for (size_t i = 0; i < std::min(a.size(), b.size()); ++i)
          if (rank(a[i]) != rank(b[i]))
            return rank(a[i]) < rank(b[i]);
        return a.size() < b.size();
      });
  }
  m_sets = m_setsByN[m_passes.front()];
  return problems;
}

std::vector<std::vector<int>> SeedFinder::buildSets(int n, std::vector<std::string>& problems) const {
  std::vector<std::vector<int>> seedLayerSets; // the result of this call (not the member m_sets)

  if (m_cfg.policy != "consecutive" && m_cfg.policy != "innermost" && m_cfg.policy != "outermost")
    problems.push_back("unknown SeedLayerPolicy '" + m_cfg.policy + "', using consecutive");

  // Sweep straight rays from the origin; the order in which a ray crosses the layers is the
  // inside-out order a track in that direction sees them.
  std::set<std::vector<int>> sets;
  const size_t len = static_cast<size_t>(n + std::max(0, m_cfg.layerSkip));
  for (double deg = 0.5; deg < 180.0; deg += 0.25) {
    std::vector<int> seq = m_geo.layersCrossedByRay(deg * M_PI / 180.0);
    if (static_cast<int>(seq.size()) < n)
      continue;
    if (m_cfg.policy == "innermost") {
      windowSubsets(seq, 0, len, n, sets);
    } else if (m_cfg.policy == "outermost") {
      std::vector<int> rev(seq.rbegin(), seq.rend());
      std::set<std::vector<int>> tmp;
      windowSubsets(rev, 0, len, n, tmp);
      for (auto s : tmp) {
        std::reverse(s.begin(), s.end());
        sets.insert(s);
      }
    } else {
      for (size_t i = 0; i + n <= seq.size(); ++i)
        windowSubsets(seq, i, len, n, sets);
    }
  }
  seedLayerSets.assign(sets.begin(), sets.end());
  return seedLayerSets;
}

// The search (best seeds first)
//   (a) seed windows inside-out; within a window seeds are extended best kink chi2 first; an
//       accepted track claims its hits immediately, later seeds using them are skipped unfitted
//   (b) hits are attached and tracks judged by the multiple-scattering-aware kink chi2
//   (c) confirmation pass: fill layers the final helix crosses without a hit, swap in better
//       hits, drop the worst kink when it is an outlier; repeat until stable
//   (d) sequential pieces of one trajectory are merged; overlapping arcs are flagged loopers
//   (e) passes: 4-hit seeds, then 3-hit seeds; the last pass claims tracks of >= MinTrackHits
namespace {
  struct TrackFit {
    double x0 = 0, y0 = 0, R = 0, circleChi2 = 0;
    Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
    std::vector<double> res;
    Helix h;
    std::vector<double> s;
    KinkResult k;
  };
} // namespace

std::vector<SECandidate> SeedFinder::run(const edm4hep::TrackerHitPlaneCollection& hits,
                                         const std::vector<std::pair<size_t, int>>& pool, SeedStats& st) const {
  st = SeedStats{};
  for (int n : m_passes)
    st.seedLayerSets += static_cast<long>(numSets(n));
  const KinkConfig kc{m_cfg.bField_T};

  std::unordered_map<size_t, int> cidOf;
  for (const auto& [g, cid] : pool)
    cidOf[g] = cid;

  // hit as the kink chi2 sees it (cached per event)
  std::unordered_map<size_t, KinkPoint> kpCache;
  auto KP = [&](size_t g) -> const KinkPoint& {
    auto it = kpCache.find(g);
    if (it != kpCache.end())
      return it->second;
    KinkPoint k;
    k.pos = posCm(hits[g]);
    k.sigT_cm = std::max(1e-4, 0.1 * hits[g].getDu());
    k.sigZ_cm = std::max(1e-4, 0.1 * hits[g].getDv());
    const auto* L = m_geo.find(cidOf.count(g) ? cidOf.at(g) : 0);
    k.xOverX0 = L ? L->xOverX0 : 0.0;
    Eigen::Vector3d o, nrm;
    if (m_plane(hits[g].getCellID(), o, nrm))
      k.normal = nrm;
    else
      k.normal =
          (L && !L->isBarrel()) ? Eigen::Vector3d::UnitZ() : Eigen::Vector3d(k.pos.x(), k.pos.y(), 0).normalized();
    return kpCache.emplace(g, k).first->second;
  };

  // full fit of an ordered hit list: circle (component's fitCircleNHits), helix, kink chi2
  auto fitTrack = [&](const std::vector<size_t>& trk, TrackFit& f) -> bool {
    if (trk.size() < 3)
      return false;
    std::vector<edm4hep::TrackerHitPlane> hv;
    std::vector<Eigen::Vector3d> pv;
    std::vector<KinkPoint> kv;
    for (size_t g : trk) {
      hv.push_back(hits[g]);
      pv.push_back(KP(g).pos);
      kv.push_back(KP(g));
    }
    if (!m_fit(hv, f.x0, f.y0, f.R, f.circleChi2, f.cov, f.res) || !(f.R > 0.0))
      return false;
    if (!helixFromHits(f.x0, f.y0, f.R, pv, f.h, &f.s))
      return false;
    f.k = kinkChi2(kv, f.R, f.h.rot, kc);
    return true;
  };
  auto maxTriplet = [](const KinkResult& k) {
    double m = 0;
    for (double v : k.perTriplet)
      m = std::max(m, v);
    return m;
  };

  std::unordered_set<size_t> claimed;
  std::map<int, LayerHits> byLayer;
  auto buildLayers = [&]() {
    byLayer.clear();
    for (const auto& [g, cid] : pool)
      if (!claimed.count(g))
        byLayer[cid].gIdx.push_back(g);
  };

  // hits on layer L near the helix crossing c, distance measured on each hit's own sensor plane
  auto nearHits = [&](const Helix& h, const Crossing& c, const LayerDescriptor& L, double W,
                      const std::vector<size_t>& trk) {
    std::vector<std::pair<double, size_t>> near;
    auto it = byLayer.find(L.compositeID);
    if (it == byLayer.end())
      return near;
    const double spread = L.isBarrel() ? (L.rMax_cm - L.rMin_cm) : (L.zMax_cm - L.zMin_cm);
    for (size_t g : it->second.gIdx) {
      if (claimed.count(g) || std::find(trk.begin(), trk.end(), g) != trk.end())
        continue;
      const Eigen::Vector3d p = KP(g).pos;
      if ((p - c.pos).norm() > W + spread)
        continue;
      double d = (p - c.pos).norm();
      Eigen::Vector3d o, nrm;
      if (m_plane(hits[g].getCellID(), o, nrm)) {
        const Crossing pc = refineToPlane(h, o, nrm, c.s);
        if (pc.ok)
          d = (p - pc.pos).norm();
      }
      if (d <= W)
        near.emplace_back(d, g);
    }
    std::sort(near.begin(), near.end());
    if (static_cast<int>(near.size()) > m_cfg.maxCandPerLayer)
      near.resize(m_cfg.maxCandPerLayer);
    return near;
  };

  // (b) extension: attach the hit with the smallest triplet kink chi2 on each next layer
  auto extendBF = [&](std::vector<size_t>& trk, int dir, std::vector<int>& visited, TrackFit& f) {
    int missed = 0;
    for (int step = 0; step < 64; ++step) {
      const double sFrom = dir > 0 ? f.s.back() : f.s.front();
      const LayerDescriptor* best = nullptr;
      Crossing bestC;
      double bestDs = std::numeric_limits<double>::infinity();
      for (const auto& L : m_geo.layers()) {
        if (std::find(visited.begin(), visited.end(), L.compositeID) != visited.end())
          continue;
        const Crossing c = crossLayer(f.h, L, sFrom, dir, m_cfg.window_cm, M_PI * f.R);
        if (!c.ok)
          continue;
        const double ds = dir * (c.s - sFrom);
        if (ds < bestDs) {
          bestDs = ds;
          bestC = c;
          best = &L;
        }
      }
      if (!best)
        break;
      visited.push_back(best->compositeID);
      const double W = m_cfg.window_cm + m_cfg.windowSlope * bestDs;
      // Judge each candidate by the chi2 increase of the REFIT that includes it, not only by
      // its kink against the current curvature: across a long gap between layers the error
      // of the current curvature alone predicts the turn wrongly by many sigma.
      std::vector<size_t> bestTrial;
      TrackFit bestFit;
      double bestD = std::numeric_limits<double>::infinity();
      const double sigR = std::sqrt(std::max(0.0, f.cov(2, 2)));
      for (const auto& [d, g] : nearHits(f.h, bestC, *best, W, trk)) {
        // first: the new hit's kink against the current curvature, including its uncertainty
        // (strict where the curvature is well measured, e.g. closely spaced layers) ...
        const size_t n = trk.size();
        const double c3 = dir > 0 ? tripletKinkChi2(KP(trk[n - 2]), KP(trk[n - 1]), KP(g), f.R, f.h.rot, kc, sigR)
                                  : tripletKinkChi2(KP(g), KP(trk[0]), KP(trk[1]), f.R, f.h.rot, kc, sigR);
        if (c3 > m_cfg.kinkMaxChi2Triplet)
          continue;
        // ... then: the chi2 increase of the refit that includes it
        std::vector<size_t> trial = trk;
        if (dir > 0)
          trial.push_back(g);
        else
          trial.insert(trial.begin(), g);
        TrackFit tf;
        if (!fitTrack(trial, tf))
          continue;
        const double dchi = tf.k.chi2 - f.k.chi2;
        if (dchi < bestD) {
          bestD = dchi;
          bestTrial = std::move(trial);
          bestFit = std::move(tf);
        }
      }
      bool accepted = false;
      if (!bestTrial.empty() && bestD <= m_cfg.kinkMaxChi2Triplet) {
        trk = std::move(bestTrial);
        f = std::move(bestFit);
        accepted = true;
      }
      if (accepted)
        missed = 0;
      else if (++missed > m_cfg.maxMissedLayers)
        break;
    }
  };

  // (c) confirmation: fill crossed-but-empty layers, swap better hits, drop outliers
  auto confirmBF = [&](std::vector<size_t>& trk, TrackFit& f) {
    for (int it = 0; it < m_cfg.confirmIterations; ++it) {
      bool changed = false;
      // fill: layers the helix crosses between the first and last hit with no hit on the track
      std::set<int> present;
      for (size_t g : trk)
        present.insert(cidOf.at(g));
      for (const auto& L : m_geo.layers()) {
        if (present.count(L.compositeID))
          continue;
        const Crossing c = crossLayer(f.h, L, f.s.front(), +1, m_cfg.window_cm, f.s.back() - f.s.front());
        if (!c.ok)
          continue;
        std::vector<size_t> bestTrial;
        TrackFit bestFit;
        double bestQ = f.k.chi2ndf() + 1.0;
        for (const auto& [d, g] : nearHits(f.h, c, L, m_cfg.window_cm, trk)) {
          const double sg = arcLengthNear(f.h, KP(g).pos, c.s);
          std::vector<size_t> trial = trk;
          size_t pos = 0;
          while (pos < f.s.size() && f.s[pos] < sg)
            ++pos;
          trial.insert(trial.begin() + pos, g);
          TrackFit tf;
          if (!fitTrack(trial, tf))
            continue;
          if (tf.k.chi2ndf() < bestQ && maxTriplet(tf.k) <= m_cfg.kinkMaxChi2Triplet) {
            bestQ = tf.k.chi2ndf();
            bestTrial = trial;
            bestFit = tf;
          }
        }
        if (!bestTrial.empty()) {
          trk = bestTrial;
          f = bestFit;
          present.insert(L.compositeID);
          ++st.confirmAdded;
          changed = true;
        }
      }
      // swap: a different hit on the same layer that lowers the chi2 clearly
      for (size_t i = 0; i < trk.size(); ++i) {
        const int cid = cidOf.at(trk[i]);
        const auto* L = m_geo.find(cid);
        auto itL = byLayer.find(cid);
        if (!L || itL == byLayer.end())
          continue;
        for (size_t g : itL->second.gIdx) {
          if (g == trk[i] || claimed.count(g) || std::find(trk.begin(), trk.end(), g) != trk.end())
            continue;
          if ((KP(g).pos - KP(trk[i]).pos).norm() > m_cfg.window_cm)
            continue;
          std::vector<size_t> trial = trk;
          trial[i] = g;
          TrackFit tf;
          if (fitTrack(trial, tf) && tf.k.chi2 < f.k.chi2 - 2.0) {
            trk = trial;
            f = tf;
            ++st.confirmSwapped;
            changed = true;
          }
        }
      }
      // outlier: remove one hit of the worst triplet, if that triplet is beyond the cut
      if (trk.size() > 3) {
        int worst = -1;
        double wv = 0;
        for (size_t i = 0; i < f.k.perTriplet.size(); ++i)
          if (f.k.perTriplet[i] > wv) {
            wv = f.k.perTriplet[i];
            worst = static_cast<int>(i);
          }
        if (worst >= 0 && wv > m_cfg.kinkMaxChi2Triplet) {
          std::vector<size_t> bestTrial;
          TrackFit bestFit;
          double bestQ = f.k.chi2ndf();
          for (int j = worst - 1; j <= worst + 1; ++j) {
            if (j < 0 || j >= static_cast<int>(trk.size()))
              continue;
            std::vector<size_t> trial = trk;
            trial.erase(trial.begin() + j);
            TrackFit tf;
            if (fitTrack(trial, tf) && tf.k.chi2ndf() < bestQ) {
              bestQ = tf.k.chi2ndf();
              bestTrial = trial;
              bestFit = tf;
            }
          }
          if (!bestTrial.empty()) {
            trk = bestTrial;
            f = bestFit;
            ++st.confirmRemoved;
            changed = true;
          }
        }
      }
      if (!changed)
        break;
    }
  };

  // every pool hit by layer, claimed or not: a hole is a crossing with no hit at all nearby
  std::unordered_map<int, std::vector<size_t>> allByLayer;
  for (const auto& [g, cid] : pool)
    allByLayer[cid].push_back(g);
  auto hitNear = [&](const LayerDescriptor& L, const Crossing& c) {
    auto it = allByLayer.find(L.compositeID);
    if (it == allByLayer.end())
      return false;
    const double spread = L.isBarrel() ? (L.rMax_cm - L.rMin_cm) : (L.zMax_cm - L.zMin_cm);
    for (size_t g : it->second)
      if ((KP(g).pos - c.pos).norm() <= m_cfg.window_cm + spread)
        return true;
    return false;
  };
  // holes: inside the track's span, then outward until it leaves the tracker (or half a turn);
  // crossings must be >= 1 cm inside the layer edges, so passing between layers is no hole
  auto countHoles = [&](const std::vector<size_t>& trk, const TrackFit& f) {
    int holes = 0;
    std::vector<int> visited;
    for (size_t g : trk)
      visited.push_back(cidOf.at(g));
    for (const auto& L : m_geo.layers()) {
      if (std::find(visited.begin(), visited.end(), L.compositeID) != visited.end())
        continue;
      const Crossing c = crossLayer(f.h, L, f.s.front(), +1, -1.0, f.s.back() - f.s.front());
      if (c.ok && !hitNear(L, c))
        ++holes;
    }
    double sFrom = f.s.back();
    for (int step = 0; step < 12 && holes <= m_cfg.shortTrackMaxHoles; ++step) {
      const LayerDescriptor* best = nullptr;
      Crossing bestC;
      double bestDs = std::numeric_limits<double>::infinity();
      for (const auto& L : m_geo.layers()) {
        if (std::find(visited.begin(), visited.end(), L.compositeID) != visited.end())
          continue;
        const Crossing c = crossLayer(f.h, L, sFrom, +1, -1.0, M_PI * f.R - (sFrom - f.s.front()));
        if (c.ok && c.s - sFrom < bestDs) {
          bestDs = c.s - sFrom;
          bestC = c;
          best = &L;
        }
      }
      if (!best)
        break;
      visited.push_back(best->compositeID);
      if (!hitNear(*best, bestC))
        ++holes;
      sFrom = bestC.s;
    }
    return holes;
  };

  std::vector<SECandidate> acc;
  std::vector<TrackFit> accFit;
  auto accept = [&](const std::vector<size_t>& trk, const TrackFit& f, int nSeed) {
    SECandidate c;
    c.gIdx = trk;
    c.chi2ndf = f.circleChi2;
    c.x0 = f.x0;
    c.y0 = f.y0;
    c.R = f.R;
    c.cov = f.cov;
    c.residuals = f.res;
    c.nSeedHits = nSeed;
    c.kinkChi2 = f.k.chi2;
    c.kinkNdf = f.k.ndf;
    c.looper = (f.s.back() - f.s.front()) > M_PI * f.R;
    for (size_t g : trk)
      claimed.insert(g);
    if (m_cfg.shadowHitWindow_cm > 0.0) {
      // unclaimed hits sitting on this trajectory (overlapping-sensor double hits) are claimed
      for (const auto& L : m_geo.layers()) {
        auto itL = allByLayer.find(L.compositeID);
        if (itL == allByLayer.end())
          continue;
        const Crossing cr = crossLayer(f.h, L, f.s.front() - 1.0, +1, 1.0, f.s.back() - f.s.front() + 2.0);
        if (!cr.ok)
          continue;
        for (size_t g : itL->second) {
          if (claimed.count(g))
            continue;
          const Eigen::Vector3d p = KP(g).pos;
          if ((p - cr.pos).norm() > 5.0)
            continue;
          Eigen::Vector3d o, nrm;
          if (!m_plane(hits[g].getCellID(), o, nrm))
            continue;
          const Crossing pc = refineToPlane(f.h, o, nrm, cr.s);
          if (pc.ok && (p - pc.pos).norm() <= m_cfg.shadowHitWindow_cm) {
            claimed.insert(g);
            ++st.shadowClaimed;
          }
        }
      }
    }
    acc.push_back(std::move(c));
    accFit.push_back(f);
  };

  // (a) + (e): passes, windows inside-out, seeds best first, immediate claiming
  for (size_t ip = 0; ip < m_passes.size(); ++ip) {
    const int n = m_passes[ip];
    const bool lastPass = (ip + 1 == m_passes.size());
    const int claimMin = lastPass ? std::max(3, m_cfg.minTrackHits) : std::max(m_cfg.claimMinHits, m_cfg.minTrackHits);
    buildLayers();
    for (const auto& set : m_setsByN.at(n)) {
      bool allPresent = true;
      for (int cid : set)
        if (!byLayer.count(cid)) {
          allPresent = false;
          break;
        }
      if (!allPresent)
        continue;

      // enumerate this window's seeds, score them by kink chi2
      struct Seed {
        std::vector<size_t> hits;
        double q;
      };
      std::vector<Seed> seeds;
      std::vector<size_t> chain;
      std::function<void(int)> grow = [&](int level) {
        if (level == n) {
          ++st.rawTriplets;
          for (size_t g : chain)
            if (claimed.count(g))
              return;
          if (n == 3 && !tripletPrecheck(KP(chain[0]).pos, KP(chain[1]).pos, KP(chain[2]).pos, m_cfg.maxZResid_cm))
            return;
          ++st.seeds;
          TrackFit f;
          if (!fitTrack(chain, f))
            return;
          for (size_t k = 1; k < f.s.size(); ++k)
            if (!(f.s[k] > f.s[k - 1]))
              return;
          if (f.k.chi2ndf() > m_cfg.seedKinkMaxChi2NDF)
            return;
          ++st.seedsAfterZ;
          seeds.push_back({chain, f.k.chi2ndf()});
          return;
        }
        const auto& L = byLayer.at(set[level]);
        if (level == 0) {
          for (size_t g : L.gIdx) {
            if (claimed.count(g))
              continue;
            chain.push_back(g);
            grow(1);
            chain.pop_back();
          }
          return;
        }
        const auto prev = KP(chain.back()).pos;
        for (size_t g : L.gIdx) {
          if (claimed.count(g))
            continue;
          if ((KP(g).pos - prev).norm() > m_cfg.maxHitDist_cm)
            continue;
          chain.push_back(g);
          if (n > 3 && chain.size() == 3 &&
              !tripletPrecheck(KP(chain[0]).pos, KP(chain[1]).pos, KP(chain[2]).pos, 2.0 * m_cfg.maxZResid_cm)) {
            chain.pop_back();
            continue;
          }
          if (n > 3 && chain.size() >= 4) {
            const size_t k = chain.size() - 1;
            if (!hitOnTrajectory(KP(chain[0]).pos, KP(chain[(k - 1) / 2]).pos, KP(chain[k - 1]).pos, KP(chain[k]).pos,
                                 m_cfg.window_cm, m_cfg.windowSlope)) {
              chain.pop_back();
              continue;
            }
          }
          grow(level + 1);
          chain.pop_back();
        }
      };
      grow(0);

      std::stable_sort(seeds.begin(), seeds.end(), [](const Seed& a, const Seed& b) { return a.q < b.q; });
      for (const auto& sd : seeds) {
        bool free = true;
        for (size_t g : sd.hits)
          if (claimed.count(g)) {
            free = false;
            break;
          }
        if (!free) {
          ++st.skippedClaimed;
          continue;
        }
        std::vector<size_t> trk = sd.hits;
        TrackFit f;
        if (!fitTrack(trk, f))
          continue;
        std::vector<int> visited(set.begin(), set.end());
        // extension direction from the seed policy: innermost seeds grow outward,
        // outermost seeds inward, consecutive seeds (anywhere) both ways
        if (m_cfg.policy != "outermost")
          extendBF(trk, +1, visited, f);
        if (m_cfg.policy != "innermost")
          extendBF(trk, -1, visited, f);
        ++st.extended;
        confirmBF(trk, f);
        if (static_cast<int>(trk.size()) < claimMin || f.k.chi2ndf() > m_cfg.trackKinkMaxChi2NDF) {
          ++st.rejectedTracks;
          continue;
        }
        if (static_cast<int>(trk.size()) <= m_cfg.shortTrackMaxHits) {
          if (m_cfg.shortTrackMinPt > 0.0 && ptFromRadius(m_cfg.bField_T, f.R) < m_cfg.shortTrackMinPt) {
            ++st.rejectedShortPt;
            continue;
          }
          if (countHoles(trk, f) > m_cfg.shortTrackMaxHoles) {
            ++st.rejectedHoles;
            continue;
          }
        }
        accept(trk, f, n);
      }
    }
  }

  // (d) merge sequential pieces of one trajectory; flag overlapping arcs as loopers
  {
    bool again = true;
    while (again) {
      again = false;
      for (size_t i = 0; i < acc.size() && !again; ++i) {
        if (acc[i].merged)
          continue;
        for (size_t j = 0; j < acc.size() && !again; ++j) {
          if (j == i || acc[j].merged || acc[j].gIdx.size() > acc[i].gIdx.size())
            continue; // i = the longer piece
          const TrackFit& fi = accFit[i];
          bool near = true, before = true, after = true;
          for (size_t g : acc[j].gIdx) {
            const Eigen::Vector3d p = KP(g).pos;
            const double sRef =
                (p - fi.h.at(fi.s.back())).norm() < (p - fi.h.at(fi.s.front())).norm() ? fi.s.back() : fi.s.front();
            const double sg = arcLengthNear(fi.h, p, sRef);
            const double ext = std::max(0.0, std::max(sg - fi.s.back(), fi.s.front() - sg));
            if ((p - fi.h.at(sg)).norm() > 3.0 * (m_cfg.window_cm + m_cfg.windowSlope * ext) || ext > M_PI * fi.R) {
              near = false;
              break;
            }
            before = before && sg < fi.s.front();
            after = after && sg > fi.s.back();
          }
          if (!near)
            continue;
          if (!before && !after) {
            acc[i].looper = acc[j].looper = true;
            continue;
          } // overlapping arcs
          std::vector<size_t> uni = before ? acc[j].gIdx : acc[i].gIdx;
          const auto& tail = before ? acc[i].gIdx : acc[j].gIdx;
          uni.insert(uni.end(), tail.begin(), tail.end());
          TrackFit fu;
          if (!fitTrack(uni, fu) || fu.k.chi2ndf() > m_cfg.trackKinkMaxChi2NDF)
            continue;
          // at the join, layers were skipped whose material the kink chi2 does not include:
          // the two triplets spanning it may exceed the single-triplet cut by up to 4x
          {
            const size_t join = (before ? acc[j].gIdx.size() : acc[i].gIdx.size()); // first hit of the second piece
            bool ok = true;
            for (size_t t = 1; t + 1 < uni.size(); ++t) {
              const bool atJoin = (t + 1 == join) || (t == join);
              if (fu.k.perTriplet[t] > (atJoin ? 4.0 : 1.0) * m_cfg.kinkMaxChi2Triplet) {
                ok = false;
                break;
              }
            }
            if (!ok)
              continue;
          }
          acc[i].gIdx = uni;
          acc[i].chi2ndf = fu.circleChi2;
          acc[i].x0 = fu.x0;
          acc[i].y0 = fu.y0;
          acc[i].R = fu.R;
          acc[i].cov = fu.cov;
          acc[i].residuals = fu.res;
          acc[i].kinkChi2 = fu.k.chi2;
          acc[i].kinkNdf = fu.k.ndf;
          acc[i].looper = acc[i].looper || (fu.s.back() - fu.s.front()) > M_PI * fu.R;
          accFit[i] = fu;
          acc[j].merged = true;
          ++st.merged;
          again = true;
        }
      }
    }
  }

  std::vector<SECandidate> out;
  for (auto& c : acc)
    if (!c.merged) {
      if (c.looper)
        ++st.loopers;
      out.push_back(std::move(c));
    }
  std::stable_sort(out.begin(), out.end(),
                   [](const SECandidate& a, const SECandidate& b) { return a.gIdx.size() > b.gIdx.size(); });
  st.candidates = static_cast<long>(out.size());
  return out;
}

} // namespace k4rt
