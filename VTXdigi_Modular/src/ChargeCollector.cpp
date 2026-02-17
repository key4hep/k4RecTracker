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
  else if (algorithm == "Debug") {
    return std::make_unique<ChargeCollector_Debug>(digitizer);
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

void ChargeCollector_SinglePixel::FillHit(const SimHitWrapper& hit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const {
  // const dd4hep::rec::Vector3D pos_global = ConvertVector(hit.simHit().getPosition());
  // const dd4hep::rec::Vector3D pos_local = GlobalToLocal(pos_global, trafoMatrix);

  // const std::array<int, 2> pix = ComputePixelIndices(pos_local, m_digitizer.m_pixelPitch, m_digitizer.m_pixelCount);

  // if (!(pix.at(0) == -1 || pix.at(1) == -1 || pix.at(0) >= static_cast<int>(m_digitizer.m_pixelCount.at(0)) || pix.at(1) >= static_cast<int>(m_digitizer.m_pixelCount.at(1))))
  //   hitMap.FillCharge(pix, hit.charge());
}

// /* -- Debug approach -- */

ChargeCollector_Debug::ChargeCollector_Debug(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_Debug constructed." << endmsg;
}

void ChargeCollector_Debug::FillHit(const SimHitWrapper& hit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const {
  const dd4hep::rec::Vector3D pos_global = ConvertVector(hit.simHit().getPosition());
  const dd4hep::rec::Vector3D pos_local = GlobalToLocal(pos_global, trafoMatrix);

  const std::pair<int, int> pix = ComputePixelIndices(pos_local, m_digitizer.m_pixelPitch, m_digitizer.m_pixelCount);
  m_digitizer.verbose() << "     - SimHit at global position       (" << pos_global.x() << ", " << pos_global.y() << ", " << pos_global.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - corresponds to local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - and pixel indices             (" << pix.first << ", " << pix.second << ")" << endmsg;

  if (pix.first == -1 || pix.second == -1 || pix.first >= static_cast<int>(m_digitizer.m_pixelCount.first) || pix.second >= static_cast<int>(m_digitizer.m_pixelCount.second)) {
    m_digitizer.warning() << "simHit local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ") is out of sensor bounds U: [" << -m_digitizer.m_sensorLength.first/2 << ", " << m_digitizer.m_sensorLength.first/2 << "], V: [" << -m_digitizer.m_sensorLength.second/2 << ", " << m_digitizer.m_sensorLength.second/2 << "]. This simHit will be skipped." << endmsg;
  }
  else {
    /* fill central pixel with 50% charge, two pixels to the right with 20% charge, pixel above with 10% charge */
    m_digitizer.verbose() << "       - Filling charge " << hit.charge() << " e." << endmsg;

    hitMap.FillCharge(pix, static_cast<int>(std::round(hit.charge() * 0.5)));
    if (pix.first + 1 < static_cast<int>(m_digitizer.m_pixelCount.first))
      hitMap.FillCharge({pix.first + 1, pix.second}, static_cast<int>(std::round(hit.charge() * 0.2)));
    if (pix.first + 2 < static_cast<int>(m_digitizer.m_pixelCount.first))
      hitMap.FillCharge({pix.first + 2, pix.second}, static_cast<int>(std::round(hit.charge() * 0.2)));
    if (pix.second + 1 < static_cast<int>(m_digitizer.m_pixelCount.second))
      hitMap.FillCharge({pix.first, pix.second + 1}, static_cast<int>(std::round(hit.charge() * 0.1)));
    
    /* split charge between this pixel and the one to the right */
    // hitMap.FillCharge(pix, static_cast<int>(std::round(hit.charge() * 0.5)));
    // if (pix.first + 1 < static_cast<int>(m_digitizer.m_pixelCount.first))
    // hitMap.FillCharge({pix.first + 1, pix.second}, static_cast<int>(std::round(hit.charge() * 0.5)));
    
    m_digitizer.verbose() << "       - Total charge collected in hitMap: " << hitMap.GetTotalCharge() << " e." << endmsg;
  }
}



// /* -- LUT approach -- */

// /* -- Drift approach -- */

} // namespace VTXdigi_tools


