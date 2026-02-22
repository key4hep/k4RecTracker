#include "../src/ChargeCollector.h"

namespace VTXdigi_tools {
using ::VTXdigi_Modular; // "unqualified name introduction from global namespace" (just so I remember what to call this in C++ speak ~ Jona)
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
  } else if (algorithm == "SinglePixel") {
    return std::make_unique<ChargeCollector_SinglePixel>(digitizer);
  }
  else if (algorithm == "Debug") {
    return std::make_unique<ChargeCollector_Debug>(digitizer);
  }
  throw std::runtime_error("Unknown ChargeCollector type: " + algorithm);
}

Path::Path(const std::shared_ptr<edm4hep::SimTrackerHit> simTrackerHit, const TGeoHMatrix& trafoMatrix, const VTXdigi_Modular& digitizer) {
  const dd4hep::rec::Vector3D simPos_global = ConvertVector(simTrackerHit->getPosition());
  simPos = GlobalToLocal(simPos_global, trafoMatrix);

  const float eps = 1e-6f; // reasonable for number O(0.01) (like sensor thickness in mm) with float precision

  double momentum_global[3] = {
    static_cast<double>(simTrackerHit->getMomentum().x), 
    static_cast<double>(simTrackerHit->getMomentum().y), 
    static_cast<double>(simTrackerHit->getMomentum().z)
  };
  double momentum_local[3];
  trafoMatrix.MasterToLocalVect(momentum_global, momentum_local);

  /* Step 1 - travel vector */
  const double scaleFactor_travel = digitizer.SensorDimensions().at(2) / std::abs(momentum_local[2]);
  travel = scaleFactor_travel * dd4hep::rec::Vector3D(momentum_local[0], momentum_local[1], momentum_local[2]); 

  /* Step 2 - entry point */
  if (abs(simPos.z()) > digitizer.SensorDimensions().at(2)/2.f + eps) {
      digitizer.warning() << "SimHit position is outside the sensor volume (local w = " << simPos.z() << " mm, sensor thickness = " << digitizer.SensorDimensions().at(2) << " mm). This should never happen. Forcing it to w=0." << endmsg;
    simPos.z() = 0.f; // ensures no divide by zero etc
  }
  float shiftDist_w;
  if (travel.z() >= 0.f)
    shiftDist_w = simPos.z() + 0.5f * digitizer.SensorDimensions().at(2);
  else
    shiftDist_w = simPos.z() - 0.5f * digitizer.SensorDimensions().at(2);
  const float scaleFactor_entry = shiftDist_w / travel.z();
  entry = simPos - scaleFactor_entry * travel;

  /* Step 3 -check that path is not much longer than the length it had in Geant4 */
  lengthG4 = simTrackerHit->getPathLength();
  if (travel.r() > 1.05f * lengthG4) {
    digitizer.debug() << "       - Shortening path length from " << static_cast<int>(travel.r()*1000) << " um to " << static_cast<int>(lengthG4*1000) << " um (the respective path length in Geant4)." << endmsg;

    /* make sure the path stays centred around the simTrackerHit position */
    const float t_simPos = ( (simPos - entry).dot(travel) ) / (travel.r() * travel.r());

    const float t_length_halved = 0.5f * lengthG4 / travel.r(); // length of the new path in terms of t [0,1] on old path, halved
    const float t_center = std::max(t_length_halved, std::min(t_simPos, 1.f - t_length_halved)); // center of new path clamped to [t_length_half, 1 - t_length_half] while not exceeding [0,1]

    const float t_min = t_center - t_length_halved;
    const float t_max = t_center + t_length_halved;

    entry = entry + t_min * travel;
    travel = (t_max - t_min) * travel;
  }

  /* Step 4 - clip path to sensor edges (in u/v) */
  std::pair<float, float> t = std::make_pair(0.f, 1.f); // parametrize path as entry + t*travel; t in [0,1]
  t = ComputePathClippingFactors(t, entry.x(), travel.x(), digitizer.SensorDimensions().at(0), digitizer);
  t = ComputePathClippingFactors(t, entry.y(), travel.y(), digitizer.SensorDimensions().at(1), digitizer);
  if (t.first != 0.f || t.second != 1.f) { 
    if (0.f <= t.first && t.first < t.second && t.second <= 1.f) {
      /* valid clipping */
      digitizer.debug() << "       - Clipping SimHitPath with t [" << t.first << ", " << t.second << "]. PathLength changed to " << static_cast<int>((t.second - t.first) * travel.r()*1000) << " um from " << static_cast<int>(travel.r()*1000) << " um" << endmsg;

      entry = entry + t.first * travel;
      travel = (t.second - t.first) * travel;
    }
    else {
      /* invalid clipping, shouldn't happen */
      digitizer.warning() << "VTXdigi_tools::Path::Path() - invalid clipping factors t = [" << t.first << ", " << t.second << "]. Path might lie completely outside the sensor." << endmsg;
      digitizer.debug() << " -> entry (" << entry.x() << ", " << entry.y() << ", " << entry.z() << ") mm, exit (" << entry.x() + travel.x() << ", " << entry.y() + travel.y() << ", " << entry.z() + travel.z() << ") mm, sensor dim. (+-" << digitizer.SensorDimensions().at(0)/2 << ", +-" << digitizer.SensorDimensions().at(1)/2 << ") mm" << endmsg;
      digitizer.debug() << " -> Path length " << static_cast<int>(travel.r()*1000) << " um, in G4 " << static_cast<int>(simTrackerHit->getPathLength()*1000) << " um" << endmsg;
    }
  }

  length = travel.r();

  digitizer.debug() << "       - Constructed path, length " << travel.r()*1000 << " um (G4-length " << lengthG4*1000 << " um), entry (" << entry.x() << ", " << entry.y() << ", " << entry.z() << ") mm, exit (" << entry.x() + travel.x() << ", " << entry.y() + travel.y() << ", " << entry.z() + travel.z() << ") mm, " << endmsg;
}

std::pair<float, float> ComputePathClippingFactors(std::pair<float,float> t, const float entry_ax, const float travel_ax, const float sensorLength_ax, const VTXdigi_Modular& digitizer) {
  /* only need the components that are parallel to the axis (u/v) that we are clipping */

  const bool positiveDir = travel_ax >= 0.f; // false -> path points in negative direction along this axis

  const float minPos = std::min(entry_ax, entry_ax + travel_ax);
  if (minPos < -0.5 * sensorLength_ax) {
    /* path extends out of sensor in negative direction*/

    const float t_clip = (-minPos - 0.5f * sensorLength_ax) / abs(travel_ax);
    if (positiveDir)
      t.first = std::max(t.first, t_clip);
    else
      t.second = std::min(t.second, 1-t_clip);

    digitizer.debug() << "       - SimHitPath extends outside sensor on NEG. side, in " << (positiveDir ? "POS" : "NEG") << ". direction. (min at " << minPos << " mm, edge at " << -0.5*sensorLength_ax << " mm) => " << -0.5*sensorLength_ax - minPos << "mm outside of the sensor, " << t_clip*100 << " percent of the path length. Clipping to [" << t.first << ", " << t.second << "]" << endmsg;
  }

  const float maxPos = std::max(entry_ax, entry_ax + travel_ax);
  if (maxPos > 0.5 * sensorLength_ax) {
    const float t_clip = (maxPos - 0.5f * sensorLength_ax) / abs(travel_ax);

    if (positiveDir)
      t.second = std::min(t.second, 1-t_clip);
    else
      t.first = std::max(t.first, t_clip);

    digitizer.debug() << "       - SimHitPath extends outside sensor on POS. side, in " << (positiveDir ? "POS" : "NEG") << ". direction. (max at " << maxPos << " mm, edge at " << 0.5*sensorLength_ax << " mm) => " << maxPos - 0.5*sensorLength_ax << "mm outside of the sensor, " << t_clip*100 << " percent of the path length. Clipping to [" << t.first << ", " << t.second << "]" << endmsg;
  }

  return t;
}



/* -- LUT approach -- */

ChargeCollector_LUT::ChargeCollector_LUT(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_LUT constructed." << endmsg;
}

void ChargeCollector_LUT::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const {

  Path path(simHit.hitPtr(), trafoMatrix, m_digitizer);

  m_digitizer.FillHistograms_fromChargeCollector_perSimHit(path.length, path.lengthG4);

  m_digitizer.debug() << "     - Filling hit map with simHit charge " << simHit.charge() << " e, path entry at (" << path.entry.x() << ", " << path.entry.y() << ", " << path.entry.z() << ") mm, travel (" << path.travel.x() << ", " << path.travel.y() << ", " << path.travel.z() << ") mm." << endmsg;
}


/* -- Single pixel approach -- */

ChargeCollector_SinglePixel::ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_SinglePixel constructed." << endmsg;
}

void ChargeCollector_SinglePixel::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const {
  const dd4hep::rec::Vector3D pos_global = ConvertVector(simHit.hitPtr()->getPosition());
  const dd4hep::rec::Vector3D pos_local = GlobalToLocal(pos_global, trafoMatrix);
  const std::pair<int, int> i_uv = ComputePixelIndices(pos_local, m_digitizer.PixelPitch(), m_digitizer.PixelCount());

  if (!(i_uv.first == -1 || i_uv.second == -1 || i_uv.first >= static_cast<int>(m_digitizer.PixelCount().first) || i_uv.second >= static_cast<int>(m_digitizer.PixelCount().second)))
    hitMap.FillCharge(i_uv, simHit.charge(), simHit.hitPtr());
}

/* -- Debug approach -- */

ChargeCollector_Debug::ChargeCollector_Debug(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_Debug constructed." << endmsg;
}

void ChargeCollector_Debug::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix) const {
  const dd4hep::rec::Vector3D pos_global = ConvertVector(simHit.hitPtr()->getPosition());
  const dd4hep::rec::Vector3D pos_local = GlobalToLocal(pos_global, trafoMatrix);
  const float charge = simHit.charge();
  const std::pair<int, int> i_uv = ComputePixelIndices(pos_local, m_digitizer.PixelPitch(), m_digitizer.PixelCount());

  m_digitizer.verbose() << "     - SimHit at global position       (" << pos_global.x() << ", " << pos_global.y() << ", " << pos_global.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - corresponds to local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - and pixel indices             (" << i_uv.first << ", " << i_uv.second << ")" << endmsg;
  if (i_uv.first == -1 || i_uv.second == -1 || i_uv.first >= static_cast<int>(m_digitizer.PixelCount().first) || i_uv.second >= static_cast<int>(m_digitizer.PixelCount().second)) {
    m_digitizer.warning() << "simHit local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ") is out of sensor bounds U: [" << -m_digitizer.SensorDimensions().at(0)/2 << ", " << m_digitizer.SensorDimensions().at(0)/2 << "], V: [" << -m_digitizer.SensorDimensions().at(1)/2 << ", " << m_digitizer.SensorDimensions().at(1)/2 << "]. This simHit will be skipped." << endmsg;
  }
  else {
    m_digitizer.verbose() << "       - Filling charge " << simHit.charge() << " e." << endmsg;

    hitMap.FillCharge(i_uv, 0.5*charge, simHit.hitPtr());
    if (i_uv.first + 1 < static_cast<int>(m_digitizer.PixelCount().first))
      hitMap.FillCharge({i_uv.first + 1, i_uv.second}, 0.2*charge, simHit.hitPtr()); 
    if (i_uv.first + 2 < static_cast<int>(m_digitizer.PixelCount().first))
      hitMap.FillCharge({i_uv.first + 2, i_uv.second}, 0.2*charge, simHit.hitPtr());
    if (i_uv.second + 1 < static_cast<int>(m_digitizer.PixelCount().second))
      hitMap.FillCharge({i_uv.first, i_uv.second + 1}, 0.1*charge, simHit.hitPtr());
      
    m_digitizer.verbose() << "       - Total charge collected in hitMap: " << hitMap.GetTotalCharge() << " e." << endmsg;
  }
}

// /* -- Drift approach -- */

} // namespace VTXdigi_tools


