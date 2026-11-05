// LayerGeometry - where each composite layer actually is.
//
// The DD4hep surfaces of each composite layer are combined into a bounding cylinder (barrel)
// or annulus (endcap disk), so that a seed can be extrapolated to "the next layer".
//
// Barrel vs endcap is decided from the surface normals; the answer of the compositeID encoding
// (below) is kept to show disagreements.
//
// Units: cm (DD4hep internal), like every other geometric quantity in the component.
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

// Layer key ("compositeID") of a hit, made by CombinatorialTrackFinder::getCompositeID:
//   barrel     100 * system + layer
//   +endcap    10000 + 100 * system + layer
//   -endcap  -(10000 + 100 * system + layer)
// The system term is there only when several subdetectors share one readout (system < 32,
// layer < 100), so a key never collides with the key of another region.
constexpr int kEndcapKeyOffset = 10000;
/// The key without its barrel/endcap part: 100 * system + layer.
inline int keyWithoutRegion(int compositeID) { return std::abs(compositeID) % kEndcapKeyOffset; }
/// The layer number inside its subdetector.
inline int layerOfKey(int compositeID) { return keyWithoutRegion(compositeID) % 100; }
/// The region from the key alone: 0 barrel, +1 +endcap, -1 -endcap.
inline int regionOfKey(int compositeID) {
  if (compositeID >= kEndcapKeyOffset)
    return 1;
  if (compositeID <= -kEndcapKeyOffset)
    return -1;
  return 0;
}

struct LayerDescriptor {
  enum class Kind { Barrel, EndcapPos, EndcapNeg };

  int compositeID = 0;
  Kind kind = Kind::Barrel;
  Kind kindFromCid = Kind::Barrel; ///< what the compositeID encoding says, for cross-checks
  double rMin_cm = 0, rMax_cm = 0, rMean_cm = 0;
  double zMin_cm = 0, zMax_cm = 0, zMean_cm = 0;
  /// Closest approach of the layer's bounding region to the origin. A single scalar that
  /// orders barrel layers and disks together from the inside out.
  double dMin_cm = 0;
  int nSurfaces = 0;     ///< sensitive surfaces only
  double xOverX0 = 0;    ///< material at normal incidence around the sensors (0 = not measured)
  int nNonSensitive = 0; ///< helper/support surfaces that were ignored

  bool isBarrel() const { return kind == Kind::Barrel; }
};

class LayerGeometry {
public:
  /// Aggregate surfaces into one descriptor per composite layer.
  void build(const std::map<int, std::vector<const dd4hep::rec::Surface*>>& surfacesByLayer);

  /// All layers, sorted inside-out by dMin_cm.
  const std::vector<LayerDescriptor>& layers() const { return m_layers; }

  /// Descriptor for a compositeID, or nullptr.
  const LayerDescriptor* find(int compositeID) const;

  /// Layers a straight ray from the origin at polar angle theta (rad) would cross, sorted
  /// by path length. Layers are treated as azimuthally complete (phi is ignored), which is
  /// right for the seeding question "which layers come first in this direction".
  std::vector<int> layersCrossedByRay(double theta) const;

  /// Human-readable problems: skipped (non-sensitive) composite IDs, normal-vs-cid
  /// disagreement, barrel radii not increasing with layer number within a system.
  /// Empty means the table looks sane.
  std::vector<std::string> validate() const;

  /// Measure each layer's material (x/X0 at normal incidence) from the DD4hep geometry: along
  /// the normal of up to nSamples of its sensors, from halfway to its inner neighbour to
  /// halfway to its outer neighbour (each side capped at maxHalfDepth_cm), so the material
  /// between layers is shared between them and nothing is counted twice.
  void measureMaterial(dd4hep::Detector& det,
                       const std::map<int, std::vector<const dd4hep::rec::Surface*>>& surfacesByLayer,
                       int nSamples = 24, double maxHalfDepth_cm = 1.0);

  bool empty() const { return m_layers.empty(); }

private:
  std::vector<LayerDescriptor> m_layers;
  std::vector<LayerDescriptor> m_skipped; ///< composite IDs with no sensitive surface
  std::map<int, size_t> m_index;          ///< compositeID -> position in m_layers
};

const char* kindName(LayerDescriptor::Kind k);

} // namespace k4rt

#endif // TRACKING_LAYERGEOMETRY_H
