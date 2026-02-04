#pragma once

#include <memory>
#include "../include/VTXdigi_Modular.h"

namespace VTXdigi_tools {

using ::VTXdigi_Modular;

/** @brief Holds position & information about path through the sensor */
struct Path;

/* -- Charge collector algorithm: LUT-based -- */

class ChargeCollector_LUT : public IChargeCollector {

public:
  explicit ChargeCollector_LUT(const VTXdigi_Modular& digitizer);
  PixelChargeMatrix Collect() const override;
};


/* -- Charge collector algorithm: Propagation-based -- */

class ChargeCollector_Drift : public IChargeCollector {

public:
  PixelChargeMatrix Collect() const override;
};


Path ComputePath(const dd4hep::rec::ISurface& surface, const edm4hep::SimTrackerHit& simHit);

} // namespace VTXdigi_tools