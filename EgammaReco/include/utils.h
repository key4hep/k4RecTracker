#ifndef EGAMMARECO_UTILS_HPP
#define EGAMMARECO_UTILS_HPP

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "DD4hep/DetElement.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Volumes.h"
#include <DDRec/DetectorData.h>

#include "edm4hep/MutableTrack.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/ClusterCollection.h"

#include <Objects/Helix.h>
#include <TDecompChol.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <marlinutil/HelixClass_double.h>

/////////////////////////////////
/// Track propagation to Ecal ///
/////////////////////////////////

/**
 * @brief Computes the track state extrapolated at the calorimeter surface.
 *
 * This function takes the projected position of the track at the calorimeter
 * (ecalProjection), the helix representation of the track at its last hit,
 * and the magnetic field along the z-axis (Bz). It extrapolates the track
 * momentum to the calorimeter position and creates a new helix at this point.
 *
 * The output is an edm4hep::TrackState filled with the helix parameters at the
 * calorimeter, including position, momentum direction, curvature, and reference point.
 *
 * @param ecalProjection CartesianVector representing the track position at the calorimeter.
 * @param helixAtLastHit HelixClass_double instance describing the track at the last hit point.
 * @param Bz Magnetic field value along the z-axis in the detector region.
 *
 * @return edm4hep::TrackState The extrapolated track state at the calorimeter surface.
 */

edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection,
                                                  const HelixClass_double& helixAtLastHit, double Bz);

/* track extrapolation to calorimeter cluster from the last track state */
edm4hep::TrackState helixAtLastTrack(const edm4hep::Track& track, double eCalBarrelInnerR, double eCalBarrelMaxZ,
                                                                double eCalEndCapInnerR, double eCalEndCapOuterR,
                                                                double eCalEndCapInnerZ, double Bz);



/* Three vector structure for convenience of calculation*/
struct ThreeVector {
    double x, y, z;

    ThreeVector() : x(0), y(0), z(0) {}
    ThreeVector(double x0, double y0, double z0) : x(x0), y(y0), z(z0) {}

    double mag() const { return std::sqrt(x*x + y*y + z*z); }
    double magXY() const { return std::sqrt(x*x + y*y); }
    double phi() const { return std::atan2(y,x); } // -pi to pi
    double theta() const { return std::acos(z / mag()); } // 0 to pi
    double eta() const { return -std::log(std::tan(theta()/2.)); }
  };


/* Calculation of phi differences */
double deltaPhi(const double p1, const double p2);

/* Computing solenoid magnetic field based on location */
double getFieldFromCompact();

/* Computing momentum from the track parameters */
ThreeVector getMomentum(const edm4hep::TrackState& state, double Bz);


/* Computing brem photon interaction point at Endcap region */
bool GetPointInZ(const double zPlane,        // z plane (endcap) coordinate
                   const ThreeVector& refPt,   // reference point
                   const ThreeVector& momVec,  // momentum
                   ThreeVector& intersectPt);

/* Computing brem photon interaction point at barrel region */
bool GetPointOnCircle(const double radius,        // (barrel) radius
                        const ThreeVector& refPt,   // reference point
                        const ThreeVector& momVec,  // momentum
                        ThreeVector& intersectPt);

/* Computing brem photon interaction point at barrel region if interaction point is at origin */
bool GetPointOnCircleForRefAtOrigin(const double radius,        // (barrel) radius
                        const ThreeVector& refPt,   // reference point
                        const ThreeVector& momVec,  // momentum
                        ThreeVector& intersectPt);


#endif // EGAMMARECO_UTILS_HPP
