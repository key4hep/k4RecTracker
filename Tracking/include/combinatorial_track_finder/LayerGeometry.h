// LayerGeometry - where each tracking layer is: the DD4hep surfaces of one layer key combined
// into a bounding cylinder (barrel) or annulus (disk), so a seed can be extrapolated to "the
// next layer". Barrel or disk comes from the surface normals (the key's answer is kept for
// cross-checks). Units: cm (DD4hep).
#ifndef TRACKING_LAYERGEOMETRY_H
#define TRACKING_LAYERGEOMETRY_H

#include <cstdlib>
#include <map>
#include <string>
#include <vector>

namespace dd4hep {
class Detector;
namespace rec {
  class Surface;
}
} // namespace dd4hep

namespace k4rt {

// Layer key of a hit (CombinatorialTrackFinder::layerKeyOf):
//   barrel 100 * system + layer, +endcap 10000 + 100 * system + layer, -endcap -(10000 + ...).
// The system term only when several subdetectors share one readout (system < 32, layer < 100),
// so keys of different regions never collide.
constexpr int kEndcapKeyOffset = 10000;
/// The key without its barrel/endcap part: 100 * system + layer.
inline int keyWithoutRegion(int layerKey) { return std::abs(layerKey) % kEndcapKeyOffset; }
/// The layer number inside its subdetector.
inline int layerNumberOfKey(int layerKey) { return keyWithoutRegion(layerKey) % 100; }
/// The region from the key alone: 0 barrel, +1 +endcap, -1 -endcap.
inline int regionFromKey(int layerKey) {
  if (layerKey >= kEndcapKeyOffset)
    return 1;
  if (layerKey <= -kEndcapKeyOffset)
    return -1;
  return 0;
}

struct LayerDescriptor {
  enum class Kind { Barrel, EndcapPos, EndcapNeg };

  int layerKey = 0;
  Kind kind = Kind::Barrel;
  Kind kindFromKey = Kind::Barrel; ///< what the layer-key encoding says, for cross-checks
  double rMin_cm = 0, rMax_cm = 0, rMean_cm = 0;
  double zMin_cm = 0, zMax_cm = 0, zMean_cm = 0;
  double closestApproachCm = 0; ///< of the bounding region to the origin: orders all layers inside-out
  int numberOfSurfaces = 0;     ///< sensitive surfaces only
  double xOverX0 = 0;           ///< material at normal incidence around the sensors (0 = not measured)
  int numberOfNonSensitive = 0; ///< helper/support surfaces that were ignored

  bool isBarrel() const { return kind == Kind::Barrel; }
  /// extent along the measured direction (cm): radial for a barrel layer, along z for a disk
  double depthCm() const { return isBarrel() ? (rMax_cm - rMin_cm) : (zMax_cm - zMin_cm); }
};

class LayerGeometry {
public:
  /// One descriptor per layer key.
  void build(const std::map<int, std::vector<const dd4hep::rec::Surface*>>& surfacesByLayer);

  /// All layers, inside-out by closestApproachCm.
  const std::vector<LayerDescriptor>& layers() const { return m_layers; }

  /// Descriptor of a layer key, or nullptr.
  const LayerDescriptor* find(int layerKey) const;

  /// Layers a straight ray from the origin at polar angle theta (rad) crosses, by path length
  /// (layers taken as complete in phi: "which layers come first in this direction").
  std::vector<int> layersCrossedByRay(double theta) const;

  /// Problems found (empty = sane): skipped layer keys, normal-vs-key disagreement, barrel radii
  /// not growing with the layer number within a system.
  std::vector<std::string> validate() const;

  /// x/X0 of each layer (normal incidence) from the DD4hep geometry, along the normal of up to
  /// numberOfSamples sensors, from halfway to its inner to halfway to its outer neighbour (each
  /// side at most maxHalfDepthCm), so material between layers is shared, not counted twice.
  void measureMaterial(dd4hep::Detector& detector,
                       const std::map<int, std::vector<const dd4hep::rec::Surface*>>& surfacesByLayer,
                       int numberOfSamples = 24, double maxHalfDepthCm = 1.0);

  bool empty() const { return m_layers.empty(); }

private:
  std::vector<LayerDescriptor> m_layers;
  std::vector<LayerDescriptor> m_skipped; ///< layer keys with no sensitive surface
  std::map<int, size_t> m_index;          ///< layer key -> position in m_layers
};

const char* kindName(LayerDescriptor::Kind kind);

} // namespace k4rt

#endif // TRACKING_LAYERGEOMETRY_H
