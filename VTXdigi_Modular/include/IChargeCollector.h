// VTXdigi_Modular/include/IChargeCollector.h
#pragma once

#include "TGeoMatrix.h"

struct VTXdigi_Modular;

namespace VTXdigi_tools {
  
  class SimHitWrapper; // forward-declare things in include/VTXdigi_tools.h
  class HitMap;

class IChargeCollector {
public:
  virtual ~IChargeCollector() = default;
  virtual void FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const = 0;

protected:
  explicit IChargeCollector(const VTXdigi_Modular& digitizer) : m_digitizer(digitizer) {}

  const VTXdigi_Modular& m_digitizer;

};

std::unique_ptr<IChargeCollector> CreateChargeCollector(const VTXdigi_Modular& digitizer, const std::string& algorithm);

} // namespace VTXdigi_tools