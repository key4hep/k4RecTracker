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

void LookupTable::SetMatrix(const Index_inPix& j_uvw, const std::vector<float>& weights) {
  if (static_cast<int>(weights.size()) != m_matrixSize*m_matrixSize)
    throw std::runtime_error("VTXdigi_tools::LookupTable::SetMatrix: weights size (" + std::to_string(weights.size()) + ") does not match matrix size (" + std::to_string(m_matrixSize*m_matrixSize) + ")");

  /* check if matrix is valid */
  float sum = 0.f;
  for (int row = 0; row < m_matrixSize; ++row) {
    for (int col = 0; col < m_matrixSize; ++col) {
      sum += weights.at(row*m_matrixSize + col);
    }
  }
  if (std::isnan(sum))
    throw std::runtime_error("VTXdigi_tools::LookupTable::SetMatrix: Charge sharing matrix for in-pixel bin (" + std::to_string(j_uvw.at(0)) + "," + std::to_string(j_uvw.at(1)) + "," + std::to_string(j_uvw.at(2)) + ") contains NaN values.");
  if (sum < 0 || sum > 1.f + 1.e-5f)
    throw std::runtime_error("VTXdigi_tools::LookupTable::SetMatrix: Charge sharing matrix for in-pixel bin (" + std::to_string(j_uvw.at(0)) + "," + std::to_string(j_uvw.at(1)) + "," + std::to_string(j_uvw.at(2)) + ") has a weight sum of " + std::to_string(sum) + ", but needs to lie in [0,1].");

  const int index = _FindIndex(j_uvw);
  for (int row = 0; row < m_matrixSize; ++row) {
    for (int col = 0; col < m_matrixSize; ++col) {
      /* weights are given in row-major order, starting at top left. 
        * We store charge sharing matrices in col-major order, starting at bottom left (lowest bin index) */
      m_matrices.at(index).at(col * m_matrixSize + row) = weights.at((m_matrixSize-1-row)*m_matrixSize + col);
    }
  }
}

void LookupTable::SetAllMatrices(const std::vector<float>& weights) {
  for (int j_u = 0; j_u < m_binCount.at(0); ++j_u) {
    for (int j_v = 0; j_v < m_binCount.at(1); ++j_v) {
      for (int j_w = 0; j_w < m_binCount.at(2); ++j_w) {
        SetMatrix({j_u, j_v, j_w}, weights);
      }
    }
  }
}

const std::vector<float>& LookupTable::GetMatrix(const Index_inPix& j_uvw) const {
  if (j_uvw.at(0) < 0 || j_uvw.at(0) >= m_binCount.at(0))
    throw std::runtime_error("VTXdigi_tools::LookupTable::GetMatrix: j_u (= " + std::to_string(j_uvw.at(0)) + ") out of range");
  if (j_uvw.at(1) < 0 || j_uvw.at(1) >= m_binCount.at(1))
    throw std::runtime_error("VTXdigi_tools::LookupTable::GetMatrix: j_v (= " + std::to_string(j_uvw.at(1)) + ") out of range");
  if (j_uvw.at(2) < 0 || j_uvw.at(2) >= m_binCount.at(2))
    throw std::runtime_error("VTXdigi_tools::LookupTable::GetMatrix: j_w (= " + std::to_string(j_uvw.at(2)) + ") out of range");
  return m_matrices[_FindIndex(j_uvw)];
}

float LookupTable::GetWeight(const Index_inPix& j_uvw, const int col, const int row) const {
  if (abs(col) > (m_matrixSize-1)/2)
    throw std::runtime_error("VTXdigi_tools::LookupTable::GetWeight: col (=" + std::to_string(col) + ") out of range of matrix size.");
  if (abs(row) > (m_matrixSize-1)/2)
    throw std::runtime_error("VTXdigi_tools::LookupTable::GetWeight: row (=" + std::to_string(row) + ") out of range of matrix size.");
  const auto& matrix = GetMatrix(j_uvw);
  return matrix[(col + (m_matrixSize-1)/2) * m_matrixSize + (row + (m_matrixSize-1)/2)];
}

int LookupTable::_FindIndex (const Index_inPix& j_uvw) const {
  if (j_uvw.at(0) < 0 || j_uvw.at(0) >= m_binCount.at(0))
    throw std::runtime_error("VTXdigi_tools::LookupTable::_FindIndex: j_u (= " + std::to_string(j_uvw.at(0)) + ") out of range");
  if (j_uvw.at(1) < 0 || j_uvw.at(1) >= m_binCount.at(1))
    throw std::runtime_error("VTXdigi_tools::LookupTable::_FindIndex: j_v (= " + std::to_string(j_uvw.at(1)) + ") out of range");
  if (j_uvw.at(2) < 0 || j_uvw.at(2) >= m_binCount.at(2))
    throw std::runtime_error("VTXdigi_tools::LookupTable::_FindIndex: j_w (= " + std::to_string(j_uvw.at(2)) + ") out of range");

  return j_uvw.at(0) + m_binCount.at(0) * (j_uvw.at(1) + m_binCount.at(1) * j_uvw.at(2)); 
}

LookupTable::LookupTable(const std::string& lutFileName, const VTXdigi_Modular& digitizer) {
  digitizer.debug() << " - Constructing LUT from file \"" << lutFileName << "\"." << endmsg;

  /* This implements parsing the default Allpix2 LUT file format.
  * See https://indico.cern.ch/event/1489052/contributions/6475539/attachments/3063712/5418424/Allpix_workshop_Lemoine.pdf (slide 10) for more info on fields in the LUT file */

  if (lutFileName.empty())
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): LUT file name is empty. A LUT file must be given to load the lookup table.");
  
  const int headerLines = 5;

  digitizer.debug() << "   - Opening LUT file \"" << lutFileName << "\"." << endmsg;
  std::ifstream lutFile(lutFileName);
  if (!lutFile.is_open())
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Could not open LUT file \"" + lutFileName + "\".");

  std::string line;
  int lineNumber = 0;

  /* Parse pixel-pitch, thickness, in-pixel bin count from header (all in 5th line) */
  for (; lineNumber < 5; ++lineNumber)
    std::getline(lutFile, line);
  std::istringstream headerStringStream(line);
  std::string headerEntry;
  std::vector<std::string> headerLineEntries;

  while (std::getline(headerStringStream, headerEntry, ' ')) {
    if (!headerEntry.empty())
      headerLineEntries.push_back(headerEntry);
  }

  if (headerLineEntries.size() != 11)
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Invalid number of entries in LUT file in 5th header line: found " + std::to_string(headerLineEntries.size()) + " entries, expected 11.");

  for (int j=0; j<3; j++) {
    m_binCount.at(j) = std::stoi(headerLineEntries.at(7+j));
  }
  digitizer.debug() << "   - found in-pixel bin count of (" << m_binCount.at(0) << ", " << m_binCount.at(1) << ", " << m_binCount.at(2) << ") from LUT file header." << endmsg;

  /* -> compare the values we just parsed to the values retrieved from the detector geometry */
  const float eps = 1e-7f; // reasonable for number O(0.01) (like sensor thickness in mm) with float precision

  const float sensorThickness = std::stof(headerLineEntries.at(0)) / 1000.f; // convert from um to mm
  if (std::abs(sensorThickness - digitizer.SensorDimensions().at(2)) > eps) 
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Sensor thickness mismatch between LUT file and detector geometry: LUT file specifies " + std::to_string(sensorThickness) + " mm, but geometry has " + std::to_string(digitizer.SensorDimensions().at(2)) + " mm.");

  const std::pair<float, float> pitch = {std::stof(headerLineEntries.at(1)) / 1000.f, std::stof(headerLineEntries.at(2)) / 1000.f};
  if (std::abs(pitch.first - digitizer.PixelPitch().first) > eps || std::abs(pitch.second - digitizer.PixelPitch().second) > eps) 
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Pixel pitch mismatch between LUT file and detector geometry: LUT file specifies (" + std::to_string(pitch.first) + ", " + std::to_string(pitch.second) + ") mm, but geometry has (" + std::to_string(digitizer.PixelPitch().first) + ", " + std::to_string(digitizer.PixelPitch().second) + ") mm.");

  digitizer.debug() << "   - Found matching pixel pitch and sensor thickness in LUT file." << endmsg;

  /* Get matrix size (5x5, 7x7, ...) from the length of the first line after the header */
  if (std::getline(lutFile, line)) {
    m_matrixSize = static_cast<int>(std::sqrt(std::count(line.begin(), line.end(), ' ') - 2)); // not very robust, but works for valid Allpix2 files. first 3 entries are bin indices
  }
  else {
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Could not read first line after header in LUT file: " + digitizer.LutFileName());
  }
  if (m_matrixSize < 3 || m_matrixSize % 2 == 0)
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Matrix size must be an odd integer >= 3, but is " + std::to_string(m_matrixSize) + ".");
  digitizer.debug() << "   - Inferred matrix size of " << m_matrixSize << " from first line." << endmsg;

  /* Set up the matrix vector */
  m_matrices.resize(m_binCount.at(0) * m_binCount.at(1) * m_binCount.at(2));
  for (auto& matrix : m_matrices) {
    matrix.resize(m_matrixSize*m_matrixSize, 0.f);
  }

  /* set up mapping from Allpix2 LUT format
  *   (row-major, starts on bottom left)
  * to the format expected by the LookupTable class 
  *   (row-major, starts on top-left) */
  std::unordered_map<int, int> indexMapping; // i: index in local format; indexMapping[i]: index in Allpix2 format
  for (int i_u = 0; i_u < m_matrixSize; i_u++) {
    for (int i_v = 0; i_v < m_matrixSize; i_v++) {
      const int i_AP2 = i_u + (m_matrixSize - 1 - i_v) * m_matrixSize;
      int i_VTXdigi = i_u + i_v * m_matrixSize;
      indexMapping[i_VTXdigi] = i_AP2;
    }
  }

  /* Now we can finally start parsing the matrix values */
  lutFile.clear();
  lutFile.seekg(0, std::ios::beg);
  for (int i=0; i<headerLines; ++i) // advance past header again
  std::getline(lutFile, line);
  
  digitizer.debug() << "   - Parsing LUT file, filling into lookup table." << endmsg;

  float matricesEntrySum = 0.f;
  lineNumber = headerLines + 1;
  while (std::getline(lutFile, line)) {
    if (line.empty() || line[0] == '#')
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Empty or comment line found in LUT file at line " + std::to_string(lineNumber+1) + ". All lines (past the 5 header lines) must contain valid matrix data.");
    
    std::istringstream stringStream(line);
    std::vector<std::string> lineEntries;
    std::string entryString;

    /* read the line & do sanity checks */
    while (std::getline(stringStream, entryString, ' ')) {
      if (!entryString.empty())
        lineEntries.push_back(entryString);
    }

    if (static_cast<int>(lineEntries.size()) != 3 + m_matrixSize*m_matrixSize)
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Invalid number of entries in LUT file at line " + std::to_string(lineNumber+1) + ": found " + std::to_string(lineEntries.size()) + " entries, but expected " + std::to_string(3 + m_matrixSize*m_matrixSize) + " (3 for bin indices, " + std::to_string(m_matrixSize*m_matrixSize) + " for matrix values).");

    /* First 3 lines are in-pixel binning indices */
    Index_inPix j_uvw({std::stoi(lineEntries[0])-1, std::stoi(lineEntries[1])-1, std::stoi(lineEntries[2])-1});// Allpix2 input is 1-indexed. Insane, I know.

    if (j_uvw.at(0) < 0 || j_uvw.at(0) >= m_binCount[0] ||
        j_uvw.at(1) < 0 || j_uvw.at(1) >= m_binCount[1] ||
        j_uvw.at(2) < 0 || j_uvw.at(2) >= m_binCount[2]) {
      throw std::runtime_error("Invalid in-pixel bin indices in LUT file at line " + std::to_string(lineNumber+1) + ": got (" + std::to_string(j_uvw.at(0)) + ", " + std::to_string(j_uvw.at(1)) + ", " + std::to_string(j_uvw.at(2)) + "), but expected ranges are [0, " + std::to_string(m_binCount[0]-1) + "], [0, " + std::to_string(m_binCount[1]-1) + "], [0, " + std::to_string(m_binCount[2]-1) + "].");
    }

    /* Parse matrix values & set it */
    std::vector<float> matrixEntries(m_matrixSize*m_matrixSize, 0.);
    float matrixEntrySum = 0.f;
    for (int i = 0; i < m_matrixSize*m_matrixSize; i++) {
      float entry = std::stof(lineEntries[3 + indexMapping[i]]); // NaN check done on sum
      matrixEntries[i] = entry;
      matrixEntrySum += entry;
    } 
    digitizer.verbose() << "   - Parsed matrix for in-pixel bin (" << j_uvw.at(0) << ", " << j_uvw.at(1) << ", " << j_uvw.at(2) << "), entry sum " << std::to_string(matrixEntrySum) << ", setting it now..." << endmsg;
    matricesEntrySum += matrixEntrySum;
    SetMatrix(j_uvw, matrixEntries);

    lineNumber++;
  } // loop over lines containing a matrix each

  if (lineNumber - (headerLines+1) != m_binCount[0] * m_binCount[1] * m_binCount[2])
    throw std::runtime_error("Invalid number of matrices loaded from file: expected " + std::to_string(m_binCount[0] * m_binCount[1] * m_binCount[2]) + " matrices (inferred from bin count in header) but found " + std::to_string(lineNumber - headerLines) + " lines.");

  matricesEntrySum /= static_cast<float>(lineNumber - headerLines);
  digitizer.info() << " - Loaded lookup table from file. Contains " << (lineNumber - headerLines) << " matrices. " << matricesEntrySum*100 << " percent of charge deposited in the sensor volume is collected by the pixels (the rest is lost, eg. due to being outside of depletion or due to trapping)." << endmsg;
}


ChargeCollector_LUT::ChargeCollector_LUT(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer), m_LUT(digitizer.LutFileName(), digitizer) {
  /* LUT is constructed in place (from file) */
  m_digitizer.info() << " - ChargeCollector_LUT constructed successfully." << endmsg;
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


