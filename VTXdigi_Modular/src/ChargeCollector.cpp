#include "../src/ChargeCollector.h"

namespace VTXdigi_tools {
using ::VTXdigi_Modular; // "unqualified name introduction from global namespace" (just so I remember what to call this in C++ speak ~ Jona)
using ::endmsg; // makes the Copilot autocomplete work better

std::unique_ptr<IChargeCollector> CreateChargeCollector(const VTXdigi_Modular& digitizer, const std::string& algorithm) {
  if (algorithm == "LookupTable") {
    // return std::make_unique<ChargeCollector_LUT>(digitizer);
    throw std::runtime_error("ChargeCollector_LUT not implemented yet.");
  } else if (algorithm == "Drift") {
    // return std::make_unique<ChargeCollector_Drift>(digitizer);
    throw std::runtime_error("ChargeCollector_Drift not implemented yet.");
  } else if (algorithm == "Fast") {
    // return std::make_unique<ChargeCollector_Fast>();
    throw std::runtime_error("ChargeCollector_Fast not implemented yet.");
  } else if (algorithm == "SinglePixel") {
    return std::make_unique<ChargeCollector_SinglePixel>(digitizer);
  }
  throw std::runtime_error("Unknown ChargeCollector type: " + algorithm);
}

struct Path  {
  dd4hep::rec::Vector3D entry;
  dd4hep::rec::Vector3D direction;
  float length;
  float lengthG4;
  int nSegments;
  float charge;
};

Path ComputePath(const dd4hep::rec::ISurface& surface, const edm4hep::SimTrackerHit& simHit) {
  Path path;
  return path;
}


// /* -- Single pixel approach -- */

ChargeCollector_SinglePixel::ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_SinglePixel constructed." << endmsg;
}

void ChargeCollector_SinglePixel::FillHit(const Hit& hit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const {
  const dd4hep::rec::Vector3D pos_global = ConvertVector(hit.simHit().getPosition());
  const dd4hep::rec::Vector3D pos_local = GlobalToLocal(pos_global, trafoMatrix);

  const std::array<int, 2> pixel_i = ComputePixelIndices(pos_local, m_digitizer.m_pixelPitch, m_digitizer.m_pixelCount);
  m_digitizer.verbose() << "     - SimHit at global position       (" << pos_global.x() << ", " << pos_global.y() << ", " << pos_global.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - corresponds to local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - and pixel indices             (" << pixel_i.at(0) << ", " << pixel_i.at(1) << ")" << endmsg;

  if (pixel_i.at(0) == -1 || pixel_i.at(1) == -1 || pixel_i.at(0) >= static_cast<int>(m_digitizer.m_pixelCount.at(0)) || pixel_i.at(1) >= static_cast<int>(m_digitizer.m_pixelCount.at(1))) {
    m_digitizer.warning() << "simHit local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ") is out of sensor bounds U: [" << -m_digitizer.m_sensorLength.at(0)/2 << ", " << m_digitizer.m_sensorLength.at(0)/2 << "], V: [" << -m_digitizer.m_sensorLength.at(1)/2 << ", " << m_digitizer.m_sensorLength.at(1)/2 << "]. This simHit will be skipped." << endmsg;
  }
  else {
    m_digitizer.verbose() << "       - Filling charge " << hit.charge() << " e." << endmsg;
    hitMap.FillCharge(pixel_i, hit.charge());
    m_digitizer.verbose() << "       - Total charge collected in hitMap: " << hitMap.GetTotalCharge() << " e." << endmsg;
  }
}


// /* -- LUT approach -- */

// /* -- Drift approach -- */

} // namespace VTXdigi_tools


