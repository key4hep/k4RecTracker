#pragma once

#include <memory>
#include "../include/VTXdigi_Modular.h"

namespace VTXdigi_tools {

class SimHitWrapper; // forward-declare things in include/VTXdigi_tools.h
class HitMap;

using ::VTXdigi_Modular;

/** @brief Holds position & information about path through the sensor */
struct Path {
  Path(const std::shared_ptr<edm4hep::SimTrackerHit> simTrackerHit, const TGeoHMatrix& trafoMatrix, const VTXdigi_Modular& digitizer);
  Path() = default;

  dd4hep::rec::Vector3D entry;
  dd4hep::rec::Vector3D travel;
  dd4hep::rec::Vector3D simPos;

  float length; // in mm
  float lengthG4;
  int nSegments;
  float charge;
};

std::pair<float, float> ComputePathClippingFactors(std::pair<float,float> t, const float entry_ax, const float travel_ax, const float sensorLength_ax, const VTXdigi_Modular& digitizer);

/* -- Charge collector algorithm: LUT-based -- */

class ChargeCollector_LUT : public IChargeCollector {
public:
  explicit ChargeCollector_LUT(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};


/* -- Charge collector algorithms for debugging -- */

class ChargeCollector_SinglePixel : public IChargeCollector {
public:
  explicit ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};

class ChargeCollector_Debug : public IChargeCollector {
public:
  explicit ChargeCollector_Debug(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};



/* -- Charge collector algorithm: Propagation-based -- */

class ChargeCollector_Drift : public IChargeCollector {
public:
  explicit ChargeCollector_Drift(const VTXdigi_Modular& digitizer);
  void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const override;
};

} // namespace VTXdigi_tools