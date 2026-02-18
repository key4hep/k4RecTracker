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

Path ComputePath(const dd4hep::rec::ISurface& surface, const edm4hep::SimTrackerHit& simTrackerHit) {


  Path path;
  return path;
}


// /* -- Single pixel approach -- */

ChargeCollector_SinglePixel::ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_SinglePixel constructed." << endmsg;
}

void ChargeCollector_SinglePixel::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const {
  const dd4hep::rec::Vector3D pos_global = ConvertVector(simHit.hitPtr()->getPosition());
  const dd4hep::rec::Vector3D pos_local = GlobalToLocal(pos_global, trafoMatrix);
  const std::pair<int, int> i_uv = ComputePixelIndices(pos_local, m_digitizer.m_pixelPitch, m_digitizer.m_pixelCount);

  if (!(i_uv.first == -1 || i_uv.second == -1 || i_uv.first >= static_cast<int>(m_digitizer.m_pixelCount.first) || i_uv.second >= static_cast<int>(m_digitizer.m_pixelCount.second)))
    hitMap.FillCharge(i_uv, simHit.charge(), simHit.hitPtr());
}

// /* -- Debug approach -- */

ChargeCollector_Debug::ChargeCollector_Debug(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_Debug constructed." << endmsg;
}

void ChargeCollector_Debug::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const {
  const dd4hep::rec::Vector3D pos_global = ConvertVector(simHit.hitPtr()->getPosition());
  const dd4hep::rec::Vector3D pos_local = GlobalToLocal(pos_global, trafoMatrix);
  const float charge = simHit.charge();
  const std::pair<int, int> i_uv = ComputePixelIndices(pos_local, m_digitizer.m_pixelPitch, m_digitizer.m_pixelCount);

  m_digitizer.verbose() << "     - SimHit at global position       (" << pos_global.x() << ", " << pos_global.y() << ", " << pos_global.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - corresponds to local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - and pixel indices             (" << i_uv.first << ", " << i_uv.second << ")" << endmsg;
  if (i_uv.first == -1 || i_uv.second == -1 || i_uv.first >= static_cast<int>(m_digitizer.m_pixelCount.first) || i_uv.second >= static_cast<int>(m_digitizer.m_pixelCount.second)) {
    m_digitizer.warning() << "simHit local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ") is out of sensor bounds U: [" << -m_digitizer.m_sensorLength.first/2 << ", " << m_digitizer.m_sensorLength.first/2 << "], V: [" << -m_digitizer.m_sensorLength.second/2 << ", " << m_digitizer.m_sensorLength.second/2 << "]. This simHit will be skipped." << endmsg;
  }
  else {
    m_digitizer.verbose() << "       - Filling charge " << simHit.charge() << " e." << endmsg;

    hitMap.FillCharge(i_uv, 0.5*charge, simHit.hitPtr());
    if (i_uv.first + 1 < static_cast<int>(m_digitizer.m_pixelCount.first))
      hitMap.FillCharge({i_uv.first + 1, i_uv.second}, 0.2*charge, simHit.hitPtr()); 
    if (i_uv.first + 2 < static_cast<int>(m_digitizer.m_pixelCount.first))
      hitMap.FillCharge({i_uv.first + 2, i_uv.second}, 0.2*charge, simHit.hitPtr());
    if (i_uv.second + 1 < static_cast<int>(m_digitizer.m_pixelCount.second))
      hitMap.FillCharge({i_uv.first, i_uv.second + 1}, 0.1*charge, simHit.hitPtr());
      
    m_digitizer.verbose() << "       - Total charge collected in hitMap: " << hitMap.GetTotalCharge() << " e." << endmsg;
  }
}



// /* -- LUT approach -- */

// /* -- Drift approach -- */

} // namespace VTXdigi_tools


