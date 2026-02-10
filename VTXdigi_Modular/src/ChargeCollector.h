#pragma once

#include <memory>
#include "../include/VTXdigi_Modular.h"

namespace VTXdigi_tools {

class Hit; // forward-declare things in include/VTXdigi_tools.h
class SensorChargeMatrix;

using ::VTXdigi_Modular;

/** @brief Holds position & information about path through the sensor */
struct Path;


/* -- Charge collector algorithm: No charge sharing, single pixel hit -- */

class ChargeCollector_SinglePixel : public IChargeCollector {
public:
  explicit ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer);
  void FillHit(const Hit& hit, SensorChargeMatrix& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};


/* -- Charge collector algorithm: LUT-based -- */

class ChargeCollector_LUT : public IChargeCollector {
public:
  explicit ChargeCollector_LUT(const VTXdigi_Modular& digitizer);
  void FillHit(const Hit& hit, SensorChargeMatrix& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};

/* -- Charge collector algorithm: Propagation-based -- */

class ChargeCollector_Drift : public IChargeCollector {
public:
  explicit ChargeCollector_Drift(const VTXdigi_Modular& digitizer);
  void FillHit(const Hit& hit, SensorChargeMatrix& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};


Path ComputePath(const dd4hep::rec::ISurface& surface, const edm4hep::SimTrackerHit& simHit);

} // namespace VTXdigi_tools