#include "LayerGeometry.h"

#include "DD4hep/Detector.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Surface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace k4rt {

const char* kindName(LayerDescriptor::Kind k) {
  switch (k) {
  case LayerDescriptor::Kind::Barrel:
    return "barrel";
  case LayerDescriptor::Kind::EndcapPos:
    return "endcap+";
  case LayerDescriptor::Kind::EndcapNeg:
    return "endcap-";
  }
  return "?";
}

namespace {
  LayerDescriptor::Kind kindFromCompositeID(int cid) {
    if (regionOfKey(cid) == 1)
      return LayerDescriptor::Kind::EndcapPos;
    if (regionOfKey(cid) == -1)
      return LayerDescriptor::Kind::EndcapNeg;
    return LayerDescriptor::Kind::Barrel;
  }
} // namespace

void LayerGeometry::build(const std::map<int, std::vector<const dd4hep::rec::Surface*>>& surfacesByLayer) {
  m_layers.clear();
  m_skipped.clear();
  m_index.clear();

  for (const auto& [cid, surfs] : surfacesByLayer) {
    LayerDescriptor d;
    d.compositeID = cid;
    d.kindFromCid = kindFromCompositeID(cid);

    double rMin = std::numeric_limits<double>::max(), rMax = 0.0;
    double zMin = std::numeric_limits<double>::max(), zMax = -std::numeric_limits<double>::max();
    double rSum = 0.0, zSum = 0.0, nzSum = 0.0;
    int nOk = 0;

    for (const auto* s : surfs) {
      if (!s)
        continue;
      // Only sensitive surfaces define a measurement layer. Helper/support surfaces carry
      // no system ID and would land in cid 0, spanning the whole tracker.
      if (!s->type().isSensitive()) {
        ++d.nNonSensitive;
        continue;
      }
      const auto o = s->origin();
      const auto u = s->u(o);
      const auto v = s->v(o);
      const auto n = s->normal(o);
      const double hu = 0.5 * s->length_along_u();
      const double hv = 0.5 * s->length_along_v();
      // Sample a 5x5 grid over the sensor, not just its corners: r^2 is convex on a plane,
      // so the point closest to the beam axis is usually INSIDE a barrel sensor, and the
      // corners alone overestimate rMin.
      for (int iu = -2; iu <= 2; ++iu) {
        for (int iv = -2; iv <= 2; ++iv) {
          const double au = 0.5 * iu * hu, av = 0.5 * iv * hv;
          const double x = o.x() + au * u.x() + av * v.x();
          const double y = o.y() + au * u.y() + av * v.y();
          const double z = o.z() + au * u.z() + av * v.z();
          const double r = std::sqrt(x * x + y * y);
          rMin = std::min(rMin, r);
          rMax = std::max(rMax, r);
          zMin = std::min(zMin, z);
          zMax = std::max(zMax, z);
        }
      }
      rSum += std::sqrt(o.x() * o.x() + o.y() * o.y());
      zSum += o.z();
      nzSum += std::abs(n.z());
      ++nOk;
    }
    if (nOk == 0) {
      // No sensitive surface at all: not a measurement layer. Remember it for the log.
      m_skipped.push_back(d);
      continue;
    }
    d.rMin_cm = rMin;
    d.rMax_cm = rMax;
    d.rMean_cm = rSum / nOk;
    d.zMin_cm = zMin;
    d.zMax_cm = zMax;
    d.zMean_cm = zSum / nOk;
    d.nSurfaces = nOk;

    // Barrel sensors face outward (normal ~ transverse), disk sensors face along z.
    const double meanAbsNz = nzSum / nOk;
    if (meanAbsNz < 0.5)
      d.kind = LayerDescriptor::Kind::Barrel;
    else
      d.kind = (d.zMean_cm >= 0.0) ? LayerDescriptor::Kind::EndcapPos : LayerDescriptor::Kind::EndcapNeg;

    // Closest approach of the bounding region to the origin.
    if (d.isBarrel()) {
      const double zc = (zMin <= 0.0 && zMax >= 0.0) ? 0.0 : std::min(std::abs(zMin), std::abs(zMax));
      d.dMin_cm = std::sqrt(rMin * rMin + zc * zc);
    } else {
      const double za = std::min(std::abs(zMin), std::abs(zMax));
      d.dMin_cm = std::sqrt(rMin * rMin + za * za);
    }
    m_layers.push_back(d);
  }

  std::stable_sort(m_layers.begin(), m_layers.end(),
                   [](const LayerDescriptor& a, const LayerDescriptor& b) { return a.dMin_cm < b.dMin_cm; });
  for (size_t i = 0; i < m_layers.size(); ++i)
    m_index[m_layers[i].compositeID] = i;
}

const LayerDescriptor* LayerGeometry::find(int compositeID) const {
  auto it = m_index.find(compositeID);
  return it == m_index.end() ? nullptr : &m_layers[it->second];
}

std::vector<int> LayerGeometry::layersCrossedByRay(double theta) const {
  const double st = std::sin(theta), ct = std::cos(theta);
  std::vector<std::pair<double, int>> hits; // (path length, cid)
  for (const auto& d : m_layers) {
    if (d.nSurfaces == 0)
      continue;
    if (d.isBarrel()) {
      if (st <= 1e-9)
        continue;
      const double path = d.rMean_cm / st;
      const double z = path * ct;
      if (z >= d.zMin_cm && z <= d.zMax_cm)
        hits.emplace_back(path, d.compositeID);
    } else {
      if (std::abs(ct) <= 1e-9)
        continue;
      if ((d.zMean_cm > 0.0) != (ct > 0.0))
        continue; // wrong hemisphere
      const double path = std::abs(d.zMean_cm / ct);
      const double r = path * st;
      if (r >= d.rMin_cm && r <= d.rMax_cm)
        hits.emplace_back(path, d.compositeID);
    }
  }
  std::sort(hits.begin(), hits.end());
  std::vector<int> out;
  out.reserve(hits.size());
  for (const auto& h : hits)
    out.push_back(h.second);
  return out;
}

std::vector<std::string> LayerGeometry::validate() const {
  std::vector<std::string> problems;
  // group barrel layers by "system" (cid / 100) to check monotonic radius in layer number
  std::map<int, std::vector<const LayerDescriptor*>> barrelBySystem;
  for (const auto& d : m_layers) {
    std::ostringstream id;
    id << "cid " << d.compositeID;
    if (d.nSurfaces == 0)
      problems.push_back(id.str() + ": no sensitive surfaces");
    if (d.kind != d.kindFromCid)
      problems.push_back(id.str() + ": normals say " + kindName(d.kind) + " but the compositeID encoding says " +
                         kindName(d.kindFromCid));
    if (d.isBarrel())
      barrelBySystem[d.compositeID / 100].push_back(&d);
  }
  for (const auto& d : m_skipped) {
    std::ostringstream m;
    m << "cid " << d.compositeID << ": skipped, " << d.nNonSensitive
      << " surfaces and none of them sensitive (helper/support surfaces)";
    problems.push_back(m.str());
  }
  for (auto& [sys, v] : barrelBySystem) {
    std::sort(v.begin(), v.end(), [](auto* a, auto* b) { return a->compositeID < b->compositeID; });
    for (size_t i = 1; i < v.size(); ++i) {
      if (v[i]->rMean_cm <= v[i - 1]->rMean_cm) {
        std::ostringstream m;
        m << "system " << sys << ": barrel cid " << v[i]->compositeID << " (r=" << v[i]->rMean_cm
          << " cm) is not outside cid " << v[i - 1]->compositeID << " (r=" << v[i - 1]->rMean_cm << " cm)";
        problems.push_back(m.str());
      }
    }
  }
  return problems;
}

void LayerGeometry::measureMaterial(dd4hep::Detector& det,
                                    const std::map<int, std::vector<const dd4hep::rec::Surface*>>& surfacesByLayer,
                                    int nSamples, double maxHalfDepth) {
  // Each layer carries the material from halfway to its inner neighbour to halfway to its
  // outer neighbour (same kind: barrel by radius with overlapping z, disks by |z| on the same
  // side with overlapping r). Material BETWEEN layers then enters the kink of the nearest hit,
  // which matters wherever it is not concentrated at the sensors (e.g. a gas detector
  // between two silicon layers). Without a neighbour on one
  // side, 1 cm is used there.
  dd4hep::rec::MaterialManager mm(det.world().volume());
  for (auto& d : m_layers) {
    const double pos = d.isBarrel() ? d.rMean_cm : std::abs(d.zMean_cm);
    double inner = -1e9, outer = 1e9;
    for (const auto& o : m_layers) {
      if (o.compositeID == d.compositeID || o.kind != d.kind)
        continue;
      bool overlap;
      double opos;
      if (d.isBarrel()) {
        overlap = o.zMin_cm < d.zMax_cm && o.zMax_cm > d.zMin_cm;
        opos = o.rMean_cm;
      } else {
        overlap = o.rMin_cm < d.rMax_cm && o.rMax_cm > d.rMin_cm;
        opos = std::abs(o.zMean_cm);
      }
      if (!overlap)
        continue;
      if (opos < pos)
        inner = std::max(inner, opos);
      else if (opos > pos)
        outer = std::min(outer, opos);
    }
    const double hIn = inner > -1e8 ? std::min(maxHalfDepth, 0.5 * (pos - inner)) : 1.0;
    const double hOut = outer < 1e8 ? std::min(maxHalfDepth, 0.5 * (outer - pos)) : 1.0;
    auto it = surfacesByLayer.find(d.compositeID);
    if (it == surfacesByLayer.end())
      continue;
    std::vector<const dd4hep::rec::Surface*> sens;
    for (const auto* s : it->second)
      if (s && s->type().isSensitive())
        sens.push_back(s);
    if (sens.empty())
      continue;
    const size_t step = std::max<size_t>(1, sens.size() / std::max(1, nSamples));
    double sum = 0.0;
    int n = 0;
    for (size_t i = 0; i < sens.size() && n < nSamples; i += step) {
      const auto o = sens[i]->origin();
      auto nv = sens[i]->normal(o);
      // orient the normal outward: radially for barrels, away from z = 0 for disks
      const double outward = d.isBarrel() ? (nv.x() * o.x() + nv.y() * o.y()) : nv.z() * o.z();
      const double sg = outward < 0 ? -1.0 : 1.0;
      const dd4hep::rec::Vector3D p0(o.x() - sg * hIn * nv.x(), o.y() - sg * hIn * nv.y(), o.z() - sg * hIn * nv.z());
      const dd4hep::rec::Vector3D p1(o.x() + sg * hOut * nv.x(), o.y() + sg * hOut * nv.y(),
                                     o.z() + sg * hOut * nv.z());
      double x = 0.0;
      for (const auto& [mat, len] : mm.materialsBetween(p0, p1)) {
        const double X0 = mat.radLength();
        if (X0 > 0.0)
          x += len / X0;
      }
      sum += x;
      ++n;
    }
    d.xOverX0 = n ? sum / n : 0.0;
  }
}

} // namespace k4rt
