#pragma once

#include <memory>
#include "DDRec/SurfaceManager.h"

struct VTXdigi_Modular;


namespace VTXdigi_tools {
  
  using ::VTXdigi_Modular; 
  
  class Hit; // forward-declare things in include/VTXdigi_tools.h
  class HitMap;

class IChargeCollector {
public:
  virtual ~IChargeCollector() = default;
  virtual void FillHit(const Hit& hit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const = 0;

protected:
  explicit IChargeCollector(const VTXdigi_Modular& digitizer) : m_digitizer(digitizer) {}

  const VTXdigi_Modular& m_digitizer;
};

std::unique_ptr<IChargeCollector> CreateChargeCollector(const VTXdigi_Modular& digitizer, const std::string& algorithm);

} // namespace VTXdigi_tools