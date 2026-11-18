#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "DD4hep/DetElement.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include <DDRec/DetectorData.h>

#include <TGeoNavigator.h>

#include "edm4hep/MutableTrack.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"

#include <TDecompChol.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <edm4hep/utils/Helix.h>

/**
 * @brief Returns the z component of the magnetic field at the origin, in Tesla.
 *
 * Returns 0 if the detector has no field defined (e.g. no geometry loaded); callers that need a field should check
 * for that.
 */
double getBzAtOrigin(const dd4hep::Detector& detector);

/////////////////////////////////
/// Track propagation to Ecal ///
/////////////////////////////////

/**
 * @brief Retrieves a unique LayeredCalorimeterData extension from the detector.
 *
 * This function queries the global DD4hep detector instance for DetElements
 * matching the given inclusion and exclusion flags. It expects exactly one
 * DetElement to be selected. If none or multiple DetElements are found, it
 * throws a std::runtime_error with an appropriate error message.
 *
 * Once the unique DetElement is found, its extension of type
 * dd4hep::rec::LayeredCalorimeterData is returned.
 *
 * @param includeFlag Bitmask flag(s) specifying which detector types to include.
 * @param excludeFlag Bitmask flag(s) specifying which detector types to exclude.
 *
 * @return Pointer to the dd4hep::rec::LayeredCalorimeterData extension of the selected DetElement.
 *
 * @throws std::runtime_error If the global detector instance is not initialized.
 * @throws std::runtime_error If the selection of DetElements is not unique or empty.
 */
dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag);

/**
 * @brief Computes the track state extrapolated at the calorimeter surface.
 *
 * @param helixAtLastHit Helix describing the track at the last hit.
 * @param length Transverse path length from the last hit to the calorimeter [mm].
 * @param Bz Magnetic field along the z-axis [T].
 * @return Track state with its reference point on the trajectory at the crossing.
 */
edm4hep::TrackState getExtrapolationAtCalorimeter(const edm4hep::utils::Helix& helixAtLastHit, double length,
                                                  double Bz);

/**
 * @brief Extrapolates the track parameters to the calorimeter surface and updates the track with this information.
 *
 * This function takes a mutable track object and magnetic field parameters to
 * extrapolate the track helix from the last measured hit position to the
 * calorimeter surface (either barrel or endcap). It calculates the track
 * momentum at the last hit, constructs a helix representation, and then
 * projects the track onto the calorimeter surfaces to find the best
 * extrapolation point.
 *
 * The extrapolation considers both the barrel and endcap calorimeter
 * geometries and selects the projection with the earliest arrival time.
 * The corresponding TrackState at the calorimeter is then attached to the
 * provided track object.
 * If neither surface has an intersection within its acceptance, no state is added.
 *
 * @param edm4hep_track MutableTrack object to be updated with calorimeter extrapolation.
 * @param Bz Magnetic field component along the z-axis (Tesla).
 * @param m_eCalBarrelInnerR Inner radius of the barrel calorimeter.
 * @param m_eCalBarrelMaxZ Maximum half-length (z) of the barrel calorimeter.
 * @param m_eCalEndCapInnerR Inner radius of the endcap calorimeter.
 * @param m_eCalEndCapOuterR Outer radius of the endcap calorimeter.
 * @param m_eCalEndCapInnerZ z-position of the inner surface of the endcap calorimeter.
 */
void FillTrackWithCalorimeterExtrapolation(edm4hep::MutableTrack& edm4hep_track, double Bz, double m_eCalBarrelInnerR,
                                           double m_eCalBarrelMaxZ, double m_eCalEndCapInnerR,
                                           double m_eCalEndCapOuterR, double m_eCalEndCapInnerZ);

/////////////////////
/// Miscellaneous ///
/////////////////////

/**
 * @brief Simple 2D point structure representing a hit in the XY plane.
 *
 * This structure is used to store the coordinates of a point in a 2D Cartesian
 * system. It is mainly intended for tracking and fitting algorithms where hits
 * are projected onto the transverse (x, y) plane.
 *
 * The structure is intentionally minimal and does not include additional
 * metadata such as errors or weights.
 */
struct Point2D_xy {
  double x;
  double y;
};

namespace ConversionUnits {
constexpr double c_mm_s = 2.998e11;       // speed of light mm/s
constexpr double a_lcio = 1e-15 * c_mm_s; // conversion constant
} // namespace ConversionUnits

/**
 * @brief RAII guard that preserves the caller's TGeoNavigator while the navigation state is modified.
 *
 * On construction, it remembers the current navigator of gGeoManager, creates a dedicated navigator and
 * makes it current. On destruction, it makes the caller's navigator current again and deletes the
 * dedicated navigator, so both the caller's navigator and the list of navigators of gGeoManager are left
 * untouched on every exit path, including exceptions. Each guard starts from a fresh navigator, so the
 * result does not depend on previous uses.
 *
 * Not thread-safe: use it only from the thread that initialized the geometry.
 *
 * A dedicated navigator is used instead of TGeoManager::DoBackupState()/DoRestoreState(), because the
 * latter only save the position in the volume hierarchy (not the current point, direction, safety, ...)
 * and TGeoNavigator has no setters for part of its state.
 */
class GeoNavigatorStateGuard {
public:
  GeoNavigatorStateGuard();
  ~GeoNavigatorStateGuard();

  GeoNavigatorStateGuard(const GeoNavigatorStateGuard&) = delete;
  GeoNavigatorStateGuard& operator=(const GeoNavigatorStateGuard&) = delete;

private:
  Int_t m_callerIndex;
  TGeoNavigator* m_fitNavigator;
};

#endif // UTILS_HPP
