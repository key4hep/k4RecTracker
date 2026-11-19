#include "utils.h"

#include "DD4hep/DD4hepUnits.h"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

double getBzAtOrigin(const dd4hep::Detector& detector) {
  const double position[3] = {0, 0, 0};
  double magneticFieldVector[3] = {0, 0, 0};
  detector.field().magneticField(position, magneticFieldVector);
  return magneticFieldVector[2] / dd4hep::tesla;
}

dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag) {

  dd4hep::rec::LayeredCalorimeterData* theExtension = 0;
  dd4hep::Detector& detector = dd4hep::Detector::getInstance();
  if (detector.detectors().empty()) {
    throw std::runtime_error("Detector is not initialized: no DetElements found.");
  }

  const std::vector<dd4hep::DetElement>& theDetectors =
      dd4hep::DetectorSelector(detector).detectors(includeFlag, excludeFlag);

  if (theDetectors.empty()) {
    // this is a standalone function, so we cannot use any message service
    std::cout << "k4RecTracker/Tracking/src/utils.cpp: No detectors found for the given selection:" << std::endl;
    std::cout << "  includeFlag: " << dd4hep::DetType(includeFlag) << " excludeFlag: " << dd4hep::DetType(excludeFlag)
              << std::endl;

    return nullptr;
  }

  int debug_lvl = 0;
  if (debug_lvl > 0) {
    std::cout << " getExtension :  includeFlag: " << dd4hep::DetType(includeFlag)
              << " excludeFlag: " << dd4hep::DetType(excludeFlag) << "  found : " << theDetectors.size()
              << "  - first det: " << theDetectors.at(0).name() << std::endl;
  }

  if (theDetectors.size() != 1) {

    std::stringstream es;
    es << " getExtension: selection is not unique - includeFlag: " << dd4hep::DetType(includeFlag)
       << " excludeFlag: " << dd4hep::DetType(excludeFlag) << " --- found detectors : ";
    for (unsigned i = 0, N = theDetectors.size(); i < N; ++i) {
      es << theDetectors.at(i).name() << ", ";
    }

    throw std::runtime_error(es.str());
  }

  theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

  return theExtension;
}

edm4hep::TrackState getExtrapolationAtCalorimeter(const edm4hep::utils::Helix& helixAtLastHit, double length,
                                                  double Bz) {
  const auto position = helixAtLastHit.positionAt(length);
  const auto momentum = helixAtLastHit.momentumAt(length);
  const edm4hep::utils::Helix helixAtCalorimeter(position, momentum, helixAtLastHit.getCharge(), Bz, position);
  return helixAtCalorimeter.toTrackState(edm4hep::TrackState::AtCalorimeter, position);
}

void FillTrackWithCalorimeterExtrapolation(edm4hep::MutableTrack& edm4hep_track, double Bz, double m_eCalBarrelInnerR,
                                           double m_eCalBarrelMaxZ, double m_eCalEndCapInnerR,
                                           double m_eCalEndCapOuterR, double m_eCalEndCapInnerZ) {
  const auto states = edm4hep_track.getTrackStates();
  const auto lastHit =
      std::ranges::find_if(states, [](const auto& state) { return state.location == edm4hep::TrackState::AtLastHit; });
  if (lastHit == states.end()) {
    return;
  }

  try {
    const auto helix = edm4hep::utils::Helix::fromTrackState(*lastHit, Bz);
    std::optional<double> bestLength;

    if (m_eCalEndCapInnerR > 0.) {
      const double zEndCap = std::copysign(m_eCalEndCapInnerZ, helix.getMomentum().z);
      if (const auto length = helix.lengthToZ(zEndCap)) {
        const auto position = helix.positionAt(*length);
        const double r = std::hypot(position.x, position.y);
        if (r >= m_eCalEndCapInnerR && r <= m_eCalEndCapOuterR) {
          bestLength = length;
        }
      }
    }

    if (m_eCalBarrelInnerR > 0.) {
      if (const auto length = helix.lengthToRadius(m_eCalBarrelInnerR);
          length && (!bestLength || *length < *bestLength)) {
        if (std::abs(helix.positionAt(*length).z) <= m_eCalBarrelMaxZ) {
          bestLength = length;
        }
      }
    }

    if (bestLength) {
      edm4hep_track.addToTrackStates(getExtrapolationAtCalorimeter(helix, *bestLength, Bz));
    }
  } catch (const std::invalid_argument&) {
    // No extrapolation can be added for an invalid or singular fitted state.
  }
}

GeoNavigatorStateGuard::GeoNavigatorStateGuard()
    : m_callerIndex(gGeoManager->GetListOfNavigators()->IndexOf(gGeoManager->GetCurrentNavigator())),
      m_fitNavigator(gGeoManager->AddNavigator()) {
  gGeoManager->SetCurrentNavigator(gGeoManager->GetListOfNavigators()->IndexOf(m_fitNavigator));
}

GeoNavigatorStateGuard::~GeoNavigatorStateGuard() {
  gGeoManager->SetCurrentNavigator(m_callerIndex);
  gGeoManager->RemoveNavigator(m_fitNavigator);
}
