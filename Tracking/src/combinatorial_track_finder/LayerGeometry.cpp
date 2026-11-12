#include "LayerGeometry.h"

#include "DD4hep/Detector.h"
#include "DDRec/MaterialManager.h"
#include "DDRec/Surface.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace k4rt {

const char* kindName(LayerDescriptor::Kind kind) {
  switch (kind) {
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
  LayerDescriptor::Kind kindFromLayerKey(int layerKey) {
    if (regionFromKey(layerKey) == 1)
      return LayerDescriptor::Kind::EndcapPos;
    if (regionFromKey(layerKey) == -1)
      return LayerDescriptor::Kind::EndcapNeg;
    return LayerDescriptor::Kind::Barrel;
  }
} // namespace

void LayerGeometry::build(const std::map<int, std::vector<const dd4hep::rec::Surface*>>& surfacesByLayer) {
  m_layers.clear();
  m_skipped.clear();
  m_index.clear();

  for (const auto& [layerKey, surfaces] : surfacesByLayer) {
    LayerDescriptor layer;
    layer.layerKey = layerKey;
    layer.kindFromKey = kindFromLayerKey(layerKey);

    double rMin = std::numeric_limits<double>::max(), rMax = 0.0;
    double zMin = std::numeric_limits<double>::max(), zMax = -std::numeric_limits<double>::max();
    double sumR = 0.0, sumZ = 0.0, sumAbsoluteNormalZ = 0.0;
    int sensitiveCount = 0;

    for (const auto* surface : surfaces) {
      if (!surface)
        continue;
      // only sensitive surfaces are measurement layers (support surfaces have no system ID and
      // would land in layer key 0, spanning the whole tracker)
      if (!surface->type().isSensitive()) {
        ++layer.numberOfNonSensitive;
        continue;
      }
      const auto origin = surface->origin();
      const auto directionU = surface->u(origin);
      const auto directionV = surface->v(origin);
      const auto normal = surface->normal(origin);
      const double halfLengthU = 0.5 * surface->length_along_u();
      const double halfLengthV = 0.5 * surface->length_along_v();
      // a 5x5 grid over the sensor, not only its corners: the point nearest the beam is usually
      // inside a barrel sensor (r^2 is convex on a plane)
      for (int stepU = -2; stepU <= 2; ++stepU) {
        for (int stepV = -2; stepV <= 2; ++stepV) {
          const double alongU = 0.5 * stepU * halfLengthU, alongV = 0.5 * stepV * halfLengthV;
          const double x = origin.x() + alongU * directionU.x() + alongV * directionV.x();
          const double y = origin.y() + alongU * directionU.y() + alongV * directionV.y();
          const double z = origin.z() + alongU * directionU.z() + alongV * directionV.z();
          const double transverseRadius = std::sqrt(x * x + y * y);
          rMin = std::min(rMin, transverseRadius);
          rMax = std::max(rMax, transverseRadius);
          zMin = std::min(zMin, z);
          zMax = std::max(zMax, z);
        }
      }
      sumR += std::sqrt(origin.x() * origin.x() + origin.y() * origin.y());
      sumZ += origin.z();
      sumAbsoluteNormalZ += std::abs(normal.z());
      ++sensitiveCount;
    }
    if (sensitiveCount == 0) {
      m_skipped.push_back(layer); // not a measurement layer, kept for the log
      continue;
    }
    layer.rMin_cm = rMin;
    layer.rMax_cm = rMax;
    layer.rMean_cm = sumR / sensitiveCount;
    layer.zMin_cm = zMin;
    layer.zMax_cm = zMax;
    layer.zMean_cm = sumZ / sensitiveCount;
    layer.numberOfSurfaces = sensitiveCount;

    // barrel sensors face outward (normal ~ transverse), disk sensors along z
    const double meanAbsoluteNormalZ = sumAbsoluteNormalZ / sensitiveCount;
    if (meanAbsoluteNormalZ < 0.5)
      layer.kind = LayerDescriptor::Kind::Barrel;
    else
      layer.kind = (layer.zMean_cm >= 0.0) ? LayerDescriptor::Kind::EndcapPos : LayerDescriptor::Kind::EndcapNeg;

    // closest approach of the bounding region to the origin
    if (layer.isBarrel()) {
      const double nearestZ = (zMin <= 0.0 && zMax >= 0.0) ? 0.0 : std::min(std::abs(zMin), std::abs(zMax));
      layer.closestApproachCm = std::sqrt(rMin * rMin + nearestZ * nearestZ);
    } else {
      const double nearestZ = std::min(std::abs(zMin), std::abs(zMax));
      layer.closestApproachCm = std::sqrt(rMin * rMin + nearestZ * nearestZ);
    }
    m_layers.push_back(layer);
  }

  std::stable_sort(m_layers.begin(), m_layers.end(), [](const LayerDescriptor& left, const LayerDescriptor& right) {
    return left.closestApproachCm < right.closestApproachCm;
  });
  for (size_t i = 0; i < m_layers.size(); ++i)
    m_index[m_layers[i].layerKey] = i;
}

const LayerDescriptor* LayerGeometry::find(int layerKey) const {
  auto entry = m_index.find(layerKey);
  return entry == m_index.end() ? nullptr : &m_layers[entry->second];
}

std::vector<int> LayerGeometry::layersCrossedByRay(double theta) const {
  const double sinTheta = std::sin(theta), cosTheta = std::cos(theta);
  std::vector<std::pair<double, int>> crossed; // (path length, layer key)
  for (const auto& layer : m_layers) {
    if (layer.numberOfSurfaces == 0)
      continue;
    if (layer.isBarrel()) {
      if (sinTheta <= 1e-9)
        continue;
      const double pathLength = layer.rMean_cm / sinTheta;
      const double z = pathLength * cosTheta;
      if (z >= layer.zMin_cm && z <= layer.zMax_cm)
        crossed.emplace_back(pathLength, layer.layerKey);
    } else {
      if (std::abs(cosTheta) <= 1e-9)
        continue;
      if ((layer.zMean_cm > 0.0) != (cosTheta > 0.0))
        continue; // wrong hemisphere
      const double pathLength = std::abs(layer.zMean_cm / cosTheta);
      const double transverseRadius = pathLength * sinTheta;
      if (transverseRadius >= layer.rMin_cm && transverseRadius <= layer.rMax_cm)
        crossed.emplace_back(pathLength, layer.layerKey);
    }
  }
  std::sort(crossed.begin(), crossed.end());
  std::vector<int> layerKeys;
  layerKeys.reserve(crossed.size());
  for (const auto& [pathLength, layerKey] : crossed)
    layerKeys.push_back(layerKey);
  return layerKeys;
}

std::vector<std::string> LayerGeometry::validate() const {
  std::vector<std::string> problems;
  // barrel layers by "system" (layer key / 100): radius must grow with the layer number
  std::map<int, std::vector<const LayerDescriptor*>> barrelBySystem;
  for (const auto& layer : m_layers) {
    std::ostringstream name;
    name << "layer key " << layer.layerKey;
    if (layer.numberOfSurfaces == 0)
      problems.push_back(name.str() + ": no sensitive surfaces");
    if (layer.kind != layer.kindFromKey)
      problems.push_back(name.str() + ": normals say " + kindName(layer.kind) + " but the layer-key encoding says " +
                         kindName(layer.kindFromKey));
    if (layer.isBarrel())
      barrelBySystem[layer.layerKey / 100].push_back(&layer);
  }
  for (const auto& layer : m_skipped) {
    std::ostringstream message;
    message << "layer key " << layer.layerKey << ": skipped, " << layer.numberOfNonSensitive
            << " surfaces and none of them sensitive (helper/support surfaces)";
    problems.push_back(message.str());
  }
  for (auto& [systemNumber, barrelLayers] : barrelBySystem) {
    std::sort(barrelLayers.begin(), barrelLayers.end(),
              [](auto* left, auto* right) { return left->layerKey < right->layerKey; });
    for (size_t i = 1; i < barrelLayers.size(); ++i) {
      if (barrelLayers[i]->rMean_cm <= barrelLayers[i - 1]->rMean_cm) {
        std::ostringstream message;
        message << "system " << systemNumber << ": barrel layer key " << barrelLayers[i]->layerKey
                << " (r=" << barrelLayers[i]->rMean_cm << " cm) is not outside layer key "
                << barrelLayers[i - 1]->layerKey << " (r=" << barrelLayers[i - 1]->rMean_cm << " cm)";
        problems.push_back(message.str());
      }
    }
  }
  return problems;
}

void LayerGeometry::measureMaterial(dd4hep::Detector& detector,
                                    const std::map<int, std::vector<const dd4hep::rec::Surface*>>& surfacesByLayer,
                                    int numberOfSamples, double maxHalfDepthCm) {
  // Each layer carries the material from halfway to its inner to halfway to its outer neighbour
  // of the same kind (barrels by radius with overlapping z, disks by |z| with overlapping r), so
  // material between layers (e.g. a gas detector between silicon layers) enters the nearest kink.
  // Without a neighbour on one side, 1 cm is used there.
  dd4hep::rec::MaterialManager materialManager(detector.world().volume());
  for (auto& layer : m_layers) {
    const double layerPosition = layer.isBarrel() ? layer.rMean_cm : std::abs(layer.zMean_cm);
    double innerNeighbour = -1e9, outerNeighbour = 1e9;
    for (const auto& other : m_layers) {
      if (other.layerKey == layer.layerKey || other.kind != layer.kind)
        continue;
      bool overlap;
      double otherPosition;
      if (layer.isBarrel()) {
        overlap = other.zMin_cm < layer.zMax_cm && other.zMax_cm > layer.zMin_cm;
        otherPosition = other.rMean_cm;
      } else {
        overlap = other.rMin_cm < layer.rMax_cm && other.rMax_cm > layer.rMin_cm;
        otherPosition = std::abs(other.zMean_cm);
      }
      if (!overlap)
        continue;
      if (otherPosition < layerPosition)
        innerNeighbour = std::max(innerNeighbour, otherPosition);
      else if (otherPosition > layerPosition)
        outerNeighbour = std::min(outerNeighbour, otherPosition);
    }
    const double depthInward =
        innerNeighbour > -1e8 ? std::min(maxHalfDepthCm, 0.5 * (layerPosition - innerNeighbour)) : 1.0;
    const double depthOutward =
        outerNeighbour < 1e8 ? std::min(maxHalfDepthCm, 0.5 * (outerNeighbour - layerPosition)) : 1.0;
    auto layerSurfaces = surfacesByLayer.find(layer.layerKey);
    if (layerSurfaces == surfacesByLayer.end())
      continue;
    std::vector<const dd4hep::rec::Surface*> sensitiveSurfaces;
    for (const auto* surface : layerSurfaces->second)
      if (surface && surface->type().isSensitive())
        sensitiveSurfaces.push_back(surface);
    if (sensitiveSurfaces.empty())
      continue;
    const size_t step = std::max<size_t>(1, sensitiveSurfaces.size() / std::max(1, numberOfSamples));
    double sum = 0.0;
    int samples = 0;
    for (size_t i = 0; i < sensitiveSurfaces.size() && samples < numberOfSamples; i += step) {
      const auto origin = sensitiveSurfaces[i]->origin();
      auto normal = sensitiveSurfaces[i]->normal(origin);
      // the normal turned outward: radially for barrels, away from z = 0 for disks
      const double outward =
          layer.isBarrel() ? (normal.x() * origin.x() + normal.y() * origin.y()) : normal.z() * origin.z();
      const double sign = outward < 0 ? -1.0 : 1.0;
      const dd4hep::rec::Vector3D start(origin.x() - sign * depthInward * normal.x(),
                                        origin.y() - sign * depthInward * normal.y(),
                                        origin.z() - sign * depthInward * normal.z());
      const dd4hep::rec::Vector3D end(origin.x() + sign * depthOutward * normal.x(),
                                      origin.y() + sign * depthOutward * normal.y(),
                                      origin.z() + sign * depthOutward * normal.z());
      double radiationLengths = 0.0;
      for (const auto& [material, length] : materialManager.materialsBetween(start, end)) {
        const double radiationLength = material.radLength();
        if (radiationLength > 0.0)
          radiationLengths += length / radiationLength;
      }
      sum += radiationLengths;
      ++samples;
    }
    layer.xOverX0 = samples ? sum / samples : 0.0;
  }
}

} // namespace k4rt
