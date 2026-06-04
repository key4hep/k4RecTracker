#include "MUONDigitizer.h"

DECLARE_COMPONENT(MUONDigitizer)

namespace {
// One readout cell that received at least one SimHit, plus the list of
// SimHits that contributed to it.
struct FiredCell {
  int                                   yIndex;  // signed strip index in y
  int                                   zIndex;  // signed strip index in z
  double                                energy;  // sum of EDep over contributing SimHits
  dd4hep::DDSegmentation::CellID        cellID;  // full cellID (used for the local transform)
  std::vector<edm4hep::SimTrackerHit>   simHits; // contributing SimHits
};
// 8-neighbour adjacency on the (y,z) strip grid.
bool areAdjacent(const FiredCell& a, const FiredCell& b) {
  const int deltaY = std::abs(a.yIndex - b.yIndex);
  const int deltaZ = std::abs(a.zIndex - b.zIndex);
  return deltaY <= 1 && deltaZ <= 1 && (deltaY + deltaZ) > 0;
}
}

MUONDigitizer::MUONDigitizer(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "MUONDigitizer") {
  declareProperty("inputSimHits", m_input_sim_hits, "Input SimTrackerHit collection");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output digitised hit collection");
  declareProperty("outputSimDigiLink", m_output_sim_digi_link,
                  "Link between digi cluster and each contributing SimHit");
}

MUONDigitizer::~MUONDigitizer() = default;

StatusCode MUONDigitizer::initialize() {
  m_randSvc = service("RndmGenSvc", false);
  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_gauss_u.initialize(m_randSvc, Rndm::Gauss(0., m_u_resolution)).isFailure() ||
      m_gauss_v.initialize(m_randSvc, Rndm::Gauss(0., m_v_resolution)).isFailure() ||
      m_flat.initialize(m_randSvc, Rndm::Flat(0., 1.)).isFailure()) {
    error() << "Couldn't initialize spatial-smearing RNGs!" << endmsg;
    return StatusCode::FAILURE;
  }
  // Optional time smearing
  if (m_t_resolution > 0.f) {
    if (m_gauss_t.initialize(m_randSvc, Rndm::Gauss(0., m_t_resolution)).isFailure()) {
      error() << "Couldn't initialize time-smearing RNG!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // check if readout exists
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder();
  m_volman  = m_geoSvc->getDetector()->volumeManager();

  //check: the readout must carry a layer field
  if (m_decoder->fieldDescription().find("layer") == std::string::npos) {
    error() << "Readout " << m_readoutName << " does not contain a 'layer' id!" << endmsg;
    return StatusCode::FAILURE;
  }

  // check surface map 
  if (m_forceHitsOntoSurface) {
    auto* theDetector = m_geoSvc->getDetector();
    auto* surfMan     = theDetector->extension<dd4hep::rec::SurfaceManager>();
    if (!surfMan) {
      error() << "No SurfaceManager extension available for the detector" << endmsg;
      return StatusCode::FAILURE;
    }
    dd4hep::DetElement det = theDetector->detector(m_detectorName);
    m_surfMap = surfMan->map(det.name());
    if (!m_surfMap) {
      error() << "Could not find surface map for detector '" << m_detectorName << "'" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  if (m_maxClusterSize < 1u) {
    error() << "maxClusterSize must be >= 1 (got " << m_maxClusterSize << ")" << endmsg;
    return StatusCode::FAILURE;
  }

  info() << "MUONDigitizer initialised:"
         << " timeWindow=" << m_timeWindow << " ns,"
         << " maxClusterSize=" << m_maxClusterSize << " strips/view (8-neighbour),"
         << " efficiency=" << m_efficiency
         << ", time smearing=" << (m_t_resolution > 0.f ? std::to_string(m_t_resolution) + " ns" : std::string("OFF"))
         << ", forceHitsOntoSurface=" << (m_forceHitsOntoSurface.value() ? "ON" : "OFF")
         << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode MUONDigitizer::finalize() { return StatusCode::SUCCESS; }

StatusCode MUONDigitizer::execute(const EventContext& ctx) const {
  const auto* sim_hits = m_input_sim_hits.get();
  auto*       digi_hits = m_output_digi_hits.createAndPut();
  auto*       links = m_output_sim_digi_link.createAndPut();
  auto*       diffX = m_simDigiDifferenceX.createAndPut();
  auto*       diffY = m_simDigiDifferenceY.createAndPut();
  auto*       diffZ = m_simDigiDifferenceZ.createAndPut();
  auto*       clusterSizesY = m_clusterSizeY.createAndPut();
  auto*       clusterSizesZ = m_clusterSizeZ.createAndPut();

  debug() << "Event " << ctx.evt() << ": " << sim_hits->size() << " input SimHits" << endmsg;

  // -------------------------------------------------------------------------
  // Step 1. Group SimHits by cellWindow = (chamberID, timeWindow).
  // "chamberID" = cellID with the y and z strip fields zeroed,
  // "timeWindow" = floor(time / m_timeWindow); SimHits in the same chamber
  // within the same 25 ns window are candidates to be merged into a cluster.
  // -------------------------------------------------------------------------
  using CellWindow = std::pair<std::uint64_t, std::int64_t>;
  std::map<CellWindow, std::vector<FiredCell>> cellWindows;

  for (const auto& sim : *sim_hits) {
    const auto cellID = sim.getCellID();
    const int  yIdx   = m_decoder->get(cellID, "y");
    const int  zIdx   = m_decoder->get(cellID, "z");

    // Build the chamber ID: same cellID, but with y and z set to 0.
    dd4hep::DDSegmentation::CellID chamberID = cellID;
    m_decoder->set(chamberID, "y", 0);
    m_decoder->set(chamberID, "z", 0);

    const auto       tWindow    = static_cast<std::int64_t>(std::floor(sim.getTime() / m_timeWindow));
    const CellWindow cellWindow{chamberID, tWindow};

    auto& cells = cellWindows[cellWindow];
    auto  it    = std::find_if(cells.begin(), cells.end(), [&](const FiredCell& c) {
      return c.yIndex == yIdx && c.zIndex == zIdx;
    });
    if (it == cells.end()) {
      cells.push_back({yIdx, zIdx, sim.getEDep(), cellID, {sim}});
    } else {
      it->energy += sim.getEDep();
      it->simHits.push_back(sim);
    }
  }

  // -------------------------------------------------------------------------
  // Step 2. For each cellWindow, run seeded 8-neighbour clustering.
  // A cluster may grow as long as the number of distinct Y-strip indices AND
  // the number of distinct Z-strip indices it covers stay <= maxClusterSize.
  // Emit one digi hit per cluster.
  // -------------------------------------------------------------------------
  for (auto& [cellWindow, cells] : cellWindows) {
    std::vector<bool> used(cells.size(), false);

    // We always pick the next seed as the highest-energy still-unused cell.
    while (true) {
      int    seedIdx = -1;
      double seedE   = -1.0;
      for (std::size_t i = 0; i < cells.size(); ++i) {
        if (!used[i] && cells[i].energy > seedE) {
          seedE   = cells[i].energy;
          seedIdx = static_cast<int>(i);
        }
      }
      if (seedIdx < 0) break;  // no unused cells left in this cellWindow

      std::vector<std::size_t> cluster = {static_cast<std::size_t>(seedIdx)};
      used[seedIdx] = true;

      // Track distinct Y- and Z-strip indices currently in the cluster.
      std::set<int> yUsed{cells[seedIdx].yIndex};
      std::set<int> zUsed{cells[seedIdx].zIndex};

      // Grow by repeatedly attaching the highest-energy fired cell that is
      // 8-neighbour-adjacent and that does not push either view beyond
      // maxClusterSize distinct strips.
      while (true) {
        int    bestIdx = -1;
        double bestE   = -1.0;
        for (std::size_t i = 0; i < cells.size(); ++i) {
          if (used[i]) continue;
          bool adjacent = false;
          for (std::size_t j : cluster) {
            if (areAdjacent(cells[i], cells[j])) { adjacent = true; break; }
          }
          if (!adjacent) continue;

          const std::size_t newY = yUsed.count(cells[i].yIndex) ? yUsed.size() : yUsed.size() + 1;
          const std::size_t newZ = zUsed.count(cells[i].zIndex) ? zUsed.size() : zUsed.size() + 1;
          if (newY > m_maxClusterSize || newZ > m_maxClusterSize) continue;

          if (cells[i].energy > bestE) {
            bestE   = cells[i].energy;
            bestIdx = static_cast<int>(i);
          }
        }
        if (bestIdx < 0) break;
        cluster.push_back(static_cast<std::size_t>(bestIdx));
        used[bestIdx] = true;
        yUsed.insert(cells[bestIdx].yIndex);
        zUsed.insert(cells[bestIdx].zIndex);
      }

      // ---- Apply per-cluster efficiency ---------------------------------
      if (m_flat.shoot() > m_efficiency) {
        continue;
      }

      // ---- Energy-weighted centroid of contributing SimHits (mm) --------
      double totalEnergy      = 0.0;
      double centroidGlobalX  = 0.0;
      double centroidGlobalY  = 0.0;
      double centroidGlobalZ  = 0.0;
      for (std::size_t idx : cluster) {
        for (const auto& simHit : cells[idx].simHits) {
          const double eDep = simHit.getEDep();
          centroidGlobalX += simHit.getPosition().x * eDep;
          centroidGlobalY += simHit.getPosition().y * eDep;
          centroidGlobalZ += simHit.getPosition().z * eDep;
          totalEnergy     += eDep;
        }
      }
      if (totalEnergy <= 0.0) {
        // Fallback: just use the seed cell's first SimHit position.
        const auto& simHit = cells[seedIdx].simHits.front();
        centroidGlobalX = simHit.getPosition().x;
        centroidGlobalY = simHit.getPosition().y;
        centroidGlobalZ = simHit.getPosition().z;
        totalEnergy     = std::max(totalEnergy, 0.0);
      } else {
        centroidGlobalX /= totalEnergy;
        centroidGlobalY /= totalEnergy;
        centroidGlobalZ /= totalEnergy;
      }

      // ---- Gaussian smear in the local frame of the seed cell -----------
      const auto& seedCellID = cells[seedIdx].cellID;
      const auto  placement  = m_volman.lookupVolumePlacement(seedCellID);
      const auto& xform      = placement.matrix();

      double simGlobalPos_cm[3] = {centroidGlobalX * dd4hep::mm,
                               centroidGlobalY * dd4hep::mm,
                               centroidGlobalZ * dd4hep::mm};
      double simLocalPos_cm[3]  = {0., 0., 0.};
      xform.MasterToLocal(simGlobalPos_cm, simLocalPos_cm);

      const double smearedLocalPos_cm[3] = {
          0.,
          simLocalPos_cm[1] + m_gauss_u.shoot() * dd4hep::mm,
          simLocalPos_cm[2] + m_gauss_v.shoot() * dd4hep::mm};

      double digiGlobalPos_cm[3] = {0., 0., 0.};
      xform.LocalToMaster(smearedLocalPos_cm, digiGlobalPos_cm);

      // ---- Optional: project smeared digi onto the actual sensitive surface
      // when it falls outside the chamber's surface bounds.
      if (m_forceHitsOntoSurface && m_surfMap) {
        auto it = m_surfMap->find(seedCellID);
        if (it != m_surfMap->end()) {
          const dd4hep::rec::ISurface* surf = it->second;
          dd4hep::rec::Vector3D smearedGlobal(digiGlobalPos_cm[0], digiGlobalPos_cm[1], digiGlobalPos_cm[2]);
          if (!surf->insideBounds(smearedGlobal)) {
            const dd4hep::rec::Vector2D lv     = surf->globalToLocal(smearedGlobal);
            const dd4hep::rec::Vector3D onSurf = surf->localToGlobal(lv);
            debug() << "  smeared hit fell outside surface bounds (distance "
                    << surf->distance(smearedGlobal) / dd4hep::mm << " mm); projecting back" << endmsg;
            digiGlobalPos_cm[0] = onSurf.x();
            digiGlobalPos_cm[1] = onSurf.y();
            digiGlobalPos_cm[2] = onSurf.z();
          }
        }
      }

      // ---- Surface u/v axes in global, from the chamber rotation --------
      const double uLocal[3] = {0., 1., 0.};
      const double vLocal[3] = {0., 0., 1.};
      double uGlobal[3], vGlobal[3];
      xform.LocalToMasterVect(uLocal, uGlobal);
      xform.LocalToMasterVect(vLocal, vGlobal);

      const float uDir[2] = {static_cast<float>(std::acos(uGlobal[2])),
                             static_cast<float>(std::atan2(uGlobal[1], uGlobal[0]))};
      const float vDir[2] = {static_cast<float>(std::acos(vGlobal[2])),
                             static_cast<float>(std::atan2(vGlobal[1], vGlobal[0]))};

      // ---- Earliest contributing SimHit time -> leading-edge hit time ---
      float hitTime = std::numeric_limits<float>::infinity();
      for (std::size_t idx : cluster) {
        for (const auto& simHit : cells[idx].simHits) {
          if (simHit.getTime() < hitTime) hitTime = simHit.getTime();
        }
      }
      if (m_t_resolution > 0.f) {
        hitTime += static_cast<float>(m_gauss_t.shoot());
      }

      // ---- Emit digi hit ------------------------------------------------
      auto digi = digi_hits->create();
      digi.setPosition({digiGlobalPos_cm[0] / dd4hep::mm,
                        digiGlobalPos_cm[1] / dd4hep::mm,
                        digiGlobalPos_cm[2] / dd4hep::mm});
      digi.setCellID(seedCellID);
      digi.setTime(hitTime);
      digi.setEDep(totalEnergy);
      digi.setU(uDir);
      digi.setV(vDir);
      digi.setDu(m_u_resolution);
      digi.setDv(m_v_resolution);

      // One link per contributing SimHit -> this digi hit.
      std::size_t nSimHits = 0;
      for (std::size_t idx : cluster) {
        for (const auto& simHit : cells[idx].simHits) {
          auto link = links->create();
          link.setFrom(digi);
          link.setTo(simHit);
          ++nSimHits;
        }
      }

      if (msgLevel(MSG::DEBUG)) {
        debug() << "  cluster: cells=" << cluster.size()
                << " (Y=" << yUsed.size() << ", Z=" << zUsed.size() << "),"
                << " SimHits=" << nSimHits << ","
                << " EDep=" << totalEnergy
                << ", pos=(" << digiGlobalPos_cm[0] / dd4hep::mm
                << ", "      << digiGlobalPos_cm[1] / dd4hep::mm
                << ", "      << digiGlobalPos_cm[2] / dd4hep::mm << ") mm,"
                << " t="     << hitTime << " ns" << endmsg;
      }

      // ---- Validation outputs ------------------------------------------
      const double residualX = (digiGlobalPos_cm[0] - simGlobalPos_cm[0]) / dd4hep::mm;
      const double residualY = (digiGlobalPos_cm[1] - simGlobalPos_cm[1]) / dd4hep::mm;
      const double residualZ = (digiGlobalPos_cm[2] - simGlobalPos_cm[2]) / dd4hep::mm;
      diffX->push_back(residualX);
      diffY->push_back(residualY);
      diffZ->push_back(residualZ);
      clusterSizesY->push_back(static_cast<int>(yUsed.size()));
      clusterSizesZ->push_back(static_cast<int>(zUsed.size()));
    }
  }

  debug() << "Event " << ctx.evt() << ": " << digi_hits->size() << " output digi clusters" << endmsg;

  // Per-event heartbeat (INFO). Emit every printFrequency events; 0 disables.
  if (m_printFrequency > 0 && (ctx.evt() % m_printFrequency == 0)) {
    info() << "Event " << ctx.evt() << ": " << sim_hits->size() << " SimHits -> "
           << digi_hits->size() << " digi clusters" << endmsg;
  }

  return StatusCode::SUCCESS;
}
