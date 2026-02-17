#pragma once

#include <memory>
#include "../include/VTXdigi_Modular.h"

namespace VTXdigi_tools {

class SimHitWrapper; // forward-declare things in include/VTXdigi_tools.h
class HitMap;

using ::VTXdigi_Modular;

/** @brief Holds position & information about path through the sensor */
struct Path;


/* -- Charge collector algorithms for debugging -- */

class ChargeCollector_SinglePixel : public IChargeCollector {
public:
  explicit ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& hit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};

class ChargeCollector_Debug : public IChargeCollector {
public:
  explicit ChargeCollector_Debug(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& hit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};


/* -- Charge collector algorithm: LUT-based -- */

class ChargeCollector_LUT : public IChargeCollector {
public:
  explicit ChargeCollector_LUT(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& hit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};

/* -- Charge collector algorithm: Propagation-based -- */

class ChargeCollector_Drift : public IChargeCollector {
public:
  explicit ChargeCollector_Drift(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& hit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};


Path ComputePath(const dd4hep::rec::ISurface& surface, const edm4hep::SimTrackerHit& simHit);

} // namespace VTXdigi_tools