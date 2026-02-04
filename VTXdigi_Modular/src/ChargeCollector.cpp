#include "../src/ChargeCollector.h"

namespace VTXdigi_tools {
using ::VTXdigi_Modular; // unqualified name introduction from global namespace (just so I remember what to call this in C++ speak ~ Jona)
using ::endmsg; // makes the Copilot autocomplete work better

std::unique_ptr<IChargeCollector> CreateChargeCollector(const VTXdigi_Modular& digitizer, const std::string& algorithm) {
  if (algorithm == "LookupTable") {
    return std::make_unique<ChargeCollector_LUT>(digitizer);
  } else if (algorithm == "Drift") {
    // return std::make_unique<ChargeCollector_Drift>(digitizer);
    throw std::runtime_error("ChargeCollector_Drift not implemented yet.");
  } else if (algorithm == "Fast") {
    // return std::make_unique<ChargeCollector_Fast>();
    throw std::runtime_error("ChargeCollector_Fast not implemented yet.");
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


/* -- LUT approach -- */

ChargeCollector_LUT::ChargeCollector_LUT(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_LUT constructed." << endmsg;
}

PixelChargeMatrix ChargeCollector_LUT::Collect() const {
  m_digitizer.verbose() << "ChargeCollector_LUT::Collect() called." << endmsg;

  PixelChargeMatrix pixelMatrix(0, 0);

  return PixelChargeMatrix(0, 0);
}


/* -- Drift approach -- */

PixelChargeMatrix ChargeCollector_Drift::Collect() const {

  return PixelChargeMatrix(0, 0);
}

} // namespace VTXdigi_tools


