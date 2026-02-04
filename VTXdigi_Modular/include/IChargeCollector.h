#pragma once

#include <memory>

struct VTXdigi_Modular;


namespace VTXdigi_tools {

class PixelChargeMatrix; // forward-declare things in include/VTXdigi_tools.h
using ::VTXdigi_Modular;

class IChargeCollector {
public:
  virtual ~IChargeCollector() = default;
  virtual PixelChargeMatrix Collect() const = 0;

protected:
  explicit IChargeCollector(const VTXdigi_Modular& digitizer) : m_digitizer(digitizer) {}

  const VTXdigi_Modular& m_digitizer;
};

std::unique_ptr<IChargeCollector> CreateChargeCollector(const VTXdigi_Modular& digitizer, const std::string& algorithm);

} // namespace VTXdigi_tools