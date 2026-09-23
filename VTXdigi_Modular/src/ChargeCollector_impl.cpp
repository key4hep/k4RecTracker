// VTXdigi_Modular/src/ChargeCollector_impl.cpp

#include "../src/ChargeCollector_impl.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

namespace VTXdigi_tools {
using ::VTXdigi_Modular; // "unqualified name introduction from global namespace" (just so I remember what to call this in C++ speak)
using ::endmsg; // makes the Copilot autocomplete work better

std::unique_ptr<IChargeCollector> CreateChargeCollector(const VTXdigi_Modular& digitizer, const std::string& algorithm) {
  std::unique_ptr<IChargeCollector> chargeCollector;

  if (algorithm == "LookupTable") {
    chargeCollector = std::make_unique<ChargeCollector_LUT>(digitizer);
  } else if (algorithm == "Drift") {
    throw std::runtime_error("ChargeCollector_Drift not implemented yet.");
  } else if (algorithm == "Fast") {
    throw std::runtime_error("ChargeCollector_Fast not implemented yet.");
  } else if (algorithm == "SinglePixel") {
    chargeCollector = std::make_unique<ChargeCollector_SinglePixel>(digitizer);
  } else if (algorithm == "Debug") {
    chargeCollector = std::make_unique<ChargeCollector_Debug>(digitizer);
  }
  else {
    throw std::runtime_error("Unknown ChargeCollector type: " + algorithm);
  }

  digitizer.info() << " - Created charge collector with algorithm " << algorithm << ". Charge collection depth center at " << chargeCollector->GetChargeCollectionDepthCenter() << " mm." << endmsg;
  return chargeCollector;
}

bool ConstructPath(Path& path, const SimHitWrapper& simHit, const TGeoHMatrix& trafoMatrix, const VTXdigi_Modular&  digitizer) {
  const double eps = 1e-12; // reasonable for number O(0.01) (like sensor thickness in mm) with float precision

  path.simPos = simHit.truthPos();

  double momentum_global[3] = {
    static_cast<double>(simHit.hitPtr()->getMomentum().x),
    static_cast<double>(simHit.hitPtr()->getMomentum().y),
    static_cast<double>(simHit.hitPtr()->getMomentum().z)
  };
  double momentum_local[3];
  trafoMatrix.MasterToLocalVect(momentum_global, momentum_local);

  /* Step 1 - travel vector */
  const double scaleFactor_travel = digitizer.ActiveVolumeDimensions().at(2) / std::abs(momentum_local[2]);
  path.travel = scaleFactor_travel * dd4hep::rec::Vector3D(momentum_local[0], momentum_local[1], momentum_local[2]);

  /* Step 2 - entry point */
  if (std::abs(path.simPos.z()) > digitizer.ActiveVolumeDimensions().at(2)/2.f + eps) {
      digitizer.warning() << "SimHit position is outside the sensor volume (local w = " << path.simPos.z() << " mm, sensor thickness = " << digitizer.ActiveVolumeDimensions().at(2) << " mm). This should never happen. Forcing it to w=0." << endmsg;
    path.simPos.z() = 0.f; // ensures no divide by zero etc
  }
  double shiftDist_w;
  if (path.travel.z() >= 0.f) {
    shiftDist_w = path.simPos.z() + 0.5f * digitizer.ActiveVolumeDimensions().at(2);
  }
  else {
    shiftDist_w = path.simPos.z() - 0.5f * digitizer.ActiveVolumeDimensions().at(2);
  }
  const double scaleFactor_entry = shiftDist_w / path.travel.z();
  path.entry = path.simPos - scaleFactor_entry * path.travel;

  /* Step 3 - clip path to sensor edges (in u/v) */
  std::array<double, 2> t = {0., 1.}; // parametrize path as entry + t*travel; t in [0,1]
  t = ComputePathClippingFactors(t, path.entry.x(), path.travel.x(), digitizer.ActiveVolumeDimensions().at(0));
  t = ComputePathClippingFactors(t, path.entry.y(), path.travel.y(), digitizer.ActiveVolumeDimensions().at(1));
  if (t[0] != 0.f || t[1] != 1.f) {
    if (0.f <= t[0] && t[0] < t[1] && t[1] <= 1.f) {
      /* valid clipping */
      digitizer.debug() << "       - Clipping SimHitPath with t [" << t[0] << ", " << t[1] << "]. PathLength changed to " << static_cast<int>((t[1] - t[0]) * path.travel.r()*1000) << " um from " << static_cast<int>(path.travel.r()*1000) << " um" << endmsg;

      path.entry = path.entry + t[0] * path.travel;
      path.travel = (t[1] - t[0]) * path.travel;
    }
    else {
      /* invalid clipping, shouldn't happen */
      digitizer.warning() << "VTXdigi_tools::Path::Path() - invalid clipping factors t = [" << t[0] << ", " << t[1] << "]. Path might lie completely outside the sensor." << endmsg;
      digitizer.debug() << " -> entry (" << path.entry.x() << ", " << path.entry.y() << ", " << path.entry.z() << ") mm, exit (" << path.entry.x() + path.travel.x() << ", " << path.entry.y() + path.travel.y() << ", " << path.entry.z() + path.travel.z() << ") mm, sensor dim. (+-" << digitizer.ActiveVolumeDimensions().at(0)/2 << ", +-" << digitizer.ActiveVolumeDimensions().at(1)/2 << ") mm" << endmsg;
      digitizer.debug() << " -> Path length " << static_cast<int>(path.travel.r()*1000) << " um, in G4 " << static_cast<int>(simHit.hitPtr()->getPathLength()*1000) << " um" << endmsg;
      return false;
    }
  }

  /* Step 4 -check that path is not much longer than the length it had in Geant4 */
  path.lengthG4 = simHit.hitPtr()->getPathLength();
  if (path.travel.r() > kPathLengthTolerance * path.lengthG4) {
    digitizer.debug() << "       - Shortening path length from " << static_cast<int>(path.travel.r()*1000) << " um to " << static_cast<int>(path.lengthG4*1000) << " um (the respective path length in Geant4)." << endmsg;

    /* make sure the path stays centred around the simTrackerHit position */
    const double t_simPos = ( (path.simPos - path.entry).dot(path.travel) ) / (path.travel.r() * path.travel.r());

    const double t_length_halved = 0.5 * path.lengthG4 / path.travel.r(); // length of the new path in terms of t [0,1] on old path, halved
    const double t_center = std::max(t_length_halved, std::min(t_simPos, 1. - t_length_halved)); // center of new path clamped to [t_length_half, 1 - t_length_half] while not exceeding [0,1]

    const double t_min = t_center - t_length_halved;
    const double t_max = t_center + t_length_halved;

    path.entry = path.entry + t_min * path.travel;
    path.travel = (t_max - t_min) * path.travel;
  }


  digitizer.debug() << "       - Constructed path, length " << path.travel.r()*1000 << " um (G4-length " << path.lengthG4*1000 << " um), entry (" << path.entry.x() << ", " << path.entry.y() << ", " << path.entry.z() << ") mm, exit (" << path.entry.x() + path.travel.x() << ", " << path.entry.y() + path.travel.y() << ", " << path.entry.z() + path.travel.z() << ") mm, " << endmsg;
  return true; // indicate valid path constructed
}

std::array<double, 2> ComputePathClippingFactors(std::array<double, 2> t, const double entry_ax, const double travel_ax, const double sensorLength_ax) {
  /* only need the components that are parallel to the axis (u/v) that we are clipping */
  const bool positiveDir = travel_ax >= 0.; // false -> path points in negative direction along this axis

  const double minPos = std::min(entry_ax, entry_ax + travel_ax);
  if (minPos < -0.5 * sensorLength_ax) {
    /* path extends out of sensor in negative direction*/

    const double t_clip = (-minPos - 0.5 * sensorLength_ax) / std::abs(travel_ax);
    if (positiveDir){
      t[0] = std::max(t[0], t_clip);
    } else {
      t[1] = std::min(t[1], 1-t_clip);
    }
  }

  const double maxPos = std::max(entry_ax, entry_ax + travel_ax);
  if (maxPos > 0.5 * sensorLength_ax) {
    const double t_clip = (maxPos - 0.5 * sensorLength_ax) / std::abs(travel_ax);

    if (positiveDir) {
      t[1] = std::min(t[1], 1-t_clip);
    } else {
      t[0] = std::max(t[0], t_clip);
    }
  }

  return t;
}



/* -- LUT approach -- */

LookupTable::LookupTable(const std::string& lutFileName, const VTXdigi_Modular& digitizer) {
  digitizer.debug() << " - Constructing LUT from file \"" << lutFileName << "\"." << endmsg;

  /* parse the LUT in Allpix Squared format
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
  if (std::abs(sensorThickness - digitizer.ActiveVolumeDimensions().at(2)) > eps) {
    if (!digitizer.LUT_ignorePitch()) {
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Sensor thickness mismatch between LUT file and detector geometry: LUT file specifies " + std::to_string(sensorThickness) + " mm, but geometry has " + std::to_string(digitizer.ActiveVolumeDimensions().at(2)) + " mm active volume thickness.");
    }
    else {
      digitizer.warning() << "Sensor thickness mismatch between LUT file and detector geometry. LUT file: " << sensorThickness << "mm, geometry (active volume thickness): " << digitizer.ActiveVolumeDimensions().at(2) << "mm. Ignored because LookupTableIgnorePitch is set to true." << endmsg;
    }
  }

  const std::array<float, 2> pitch = {std::stof(headerLineEntries.at(1)) / 1000.f, std::stof(headerLineEntries.at(2)) / 1000.f};
  if (std::abs(pitch[0] - digitizer.PixelPitch().at(0)) > eps || std::abs(pitch[1] - digitizer.PixelPitch().at(1)) > eps) {
    if (!digitizer.LUT_ignorePitch())
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Pixel pitch mismatch between LUT file and detector geometry: LUT file specifies (" + std::to_string(pitch[0]) + ", " + std::to_string(pitch[1]) + ") mm, but geometry has (" + std::to_string(digitizer.PixelPitch().at(0)) + ", " + std::to_string(digitizer.PixelPitch().at(1)) + ") mm.");
    else
      digitizer.warning() << "Pixel pitch mismatch between LUT file and detector geometry. LUT file: " << pitch[0] << "mm, geometry: " << digitizer.PixelPitch().at(0) << "mm. Ignored because LookupTableIgnorePitch is set to true." << endmsg;
  }

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
  m_matrixSize_half = (m_matrixSize - 1) / 2;
  digitizer.debug() << "   - Inferred matrix size of " << m_matrixSize << " from first line." << endmsg;

  /* Set up the matrix vector */
  m_matrices.resize(m_binCount.at(0) * m_binCount.at(1) * m_binCount.at(2) * m_matrixSize * m_matrixSize, 0.f);

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

  std::vector<float> matricesEntrySum_perWBin;
  matricesEntrySum_perWBin.resize(m_binCount.at(2), 0.f);
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

    /* First 3 entries are in-pixel binning indices */
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
      if (entry < kLutEntryMinimum)
        entry = 0.f; // avoid very small entries for performace
      matrixEntries[i] = entry;
      matrixEntrySum += entry;
    }
    digitizer.verbose() << "   - Parsed matrix for in-pixel bin (" << j_uvw.at(0) << ", " << j_uvw.at(1) << ", " << j_uvw.at(2) << "), entry sum " << std::to_string(matrixEntrySum) << ", setting it now..." << endmsg;
    if (std::isnan(matrixEntrySum))
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Charge sharing matrix for in-pixel bin (" + std::to_string(j_uvw.at(0)) + "," + std::to_string(j_uvw.at(1)) + "," + std::to_string(j_uvw.at(2)) + ") contains NaN values (sum of entries is NaN).");

    matricesEntrySum += matrixEntrySum;
    matricesEntrySum_perWBin[j_uvw.at(2)] += matrixEntrySum;
    SetMatrix(j_uvw, matrixEntries);

    lineNumber++;
  } // loop over lines containing a matrix each

  if (lineNumber - (headerLines+1) != m_binCount[0] * m_binCount[1] * m_binCount[2])
    throw std::runtime_error("Invalid number of matrices loaded from file: expected " + std::to_string(m_binCount[0] * m_binCount[1] * m_binCount[2]) + " matrices (inferred from bin count in header) but found " + std::to_string(lineNumber - headerLines) + " lines.");

  const float matrixNumber = static_cast<float>(lineNumber - headerLines);

  // From which w-level are charges collected?
  float collectedFromW = 0.f;
  for (int j_w = 0; j_w < m_binCount.at(2); j_w++) {
    const float binHeight = sensorThickness / m_binCount.at(2);
    const float w_bin_center = j_w*binHeight + 0.5f*binHeight - 0.5f*sensorThickness;
    collectedFromW += matricesEntrySum_perWBin.at(j_w) * w_bin_center;
  }
  collectedFromW /= matricesEntrySum;

  if (std::abs(collectedFromW) > 0.5f * sensorThickness) throw std::runtime_error("Invalid charge collection depth inferred from LUT file: " + std::to_string(collectedFromW) + " mm. This is outside the sensor volume, which extends from " + std::to_string(-0.5f*sensorThickness) + " mm to " + std::to_string(0.5f*sensorThickness) + " mm.");
  if (std::abs(collectedFromW) >= 1.e-5f) {
    m_chargeCollectionDepthCenter = collectedFromW; // the member is initalised to 0.f
  }

  digitizer.info() << " - Loaded lookup table from file. Matrices parsed: " << (lineNumber - headerLines) << ". Charge collected from sensitive volume: " << matricesEntrySum/matrixNumber*100 << " percent (rest is lost, eg recombination). The center of the collection region is at " << m_chargeCollectionDepthCenter << endmsg;
}

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

  for (int row = 0; row < m_matrixSize; ++row) {
    for (int col = 0; col < m_matrixSize; ++col) {
      /* weights are given in row-major order, starting at top left.
        * We store charge sharing matrices in col-major order, starting at bottom left (lowest bin index) */
      m_matrices.at(FindIndex(j_uvw, col, row)) = weights.at((m_matrixSize-1-row)*m_matrixSize + col);
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

int LookupTable::FindIndex (const Index_inPix& j, const int col, const int row) const {
  if (j[0] < 0 || j[0] >= m_binCount[0]
    || j[1] < 0 || j[1] >= m_binCount[1]
    || j[2] < 0 || j[2] >= m_binCount[2] ) {
    throw std::runtime_error("VTXdigi_tools::LookupTable::FindIndex: in-pix bin out of range");
  }
  if (col < 0 || col >= m_matrixSize || row < 0 || row >= m_matrixSize) {
    throw std::runtime_error("VTXdigi_tools::LookupTable::FindIndex: col or row out of range");
  }

  int index_matrix = j[0] + m_binCount[0] * (j[1] + m_binCount[1] * j[2]);
  int index_element = col * m_matrixSize + row;
  return index_matrix * m_matrixSize * m_matrixSize + index_element;
}



ChargeCollector_LUT::ChargeCollector_LUT(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer),
  m_LUT(digitizer.LutFileName(), digitizer),
  m_shiftTruthPos(digitizer.LUT_shiftTruthPos()) {

  /* LUT is constructed in place (from file) */
  m_chargeCollectionDepthCenter = m_LUT.GetChargeCollectionDepthCenter();

  /* fine grid for the voxel traversal, fixed per job (all sensors share pitch & dimensions, checked at init) */
  const Index_inPix binCount = m_LUT.GetBinCount();
  const std::array<float, 2> pitch = digitizer.PixelPitch();
  const std::array<size_t, 2> pixelCount = digitizer.PixelCount();
  const float activeThickness = digitizer.ActiveVolumeDimensions().at(2);

  m_cellSize = {
    pitch[0] / binCount[0],
    pitch[1] / binCount[1],
    activeThickness / binCount[2]};
  m_gridOrigin = {
    -0.5f * pitch[0] * static_cast<float>(pixelCount[0]),
    -0.5f * pitch[1] * static_cast<float>(pixelCount[1]),
    -0.5f * activeThickness};
  m_gridBinCount = {
    static_cast<int>(pixelCount[0]) * binCount[0],
    static_cast<int>(pixelCount[1]) * binCount[1],
    binCount[2]};

  /* Load charge deposition sampling parameters */
  m_meanDepositionsPerUm = digitizer.MeanDepositionsPerUm();

  std::unique_ptr<TFile> file(TFile::Open(digitizer.DepositionChargeHistogramFilename().c_str(), "READ"));
  if (!file || file->IsZombie())
    throw std::runtime_error("Could not open deposition charge histogram file " + digitizer.DepositionChargeHistogramFilename() + ", cannot continue.");

  TH1D* hist_chargeDep = file->Get<TH1D>("deposition_charge");
  if (!hist_chargeDep)
    throw std::runtime_error("Could not find histogram \"deposition_charge\" in file "+ digitizer.DepositionChargeHistogramFilename() + ", cannot continue.");

  // move ownership from TFile to this class
  hist_chargeDep->SetDirectory(nullptr);
  hist_chargeDep->ComputeIntegral(); // so that we don't compute the integral later down the line (which would be slow & I am not sure about thread safety)
  m_chargeSamplingHist.reset(hist_chargeDep);

  m_digitizer.info() << " - ChargeCollector_LUT constructed successfully." << endmsg;
}


void ChargeCollector_LUT::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix, TRandom3& randomGen) const {

  Path path;
  if (!ConstructPath(path, simHit, trafoMatrix, m_digitizer))
    return;

  if (m_shiftTruthPos) {
    MoveTruthPosition(simHit, path); // shifts the sim hit position to the depth in the sensor where most charge is collected, to get usesful residual plots.
  }

  float probablity = 0.8;
  float mean = path.travel.r() * m_meanDepositionsPerUm * 1000.f; // convert from um to mm
  int NDepositions = static_cast<int>(randomGen.Binomial(std::round(mean / probablity), probablity)); // convert from um to mm

  std::vector<float> depositionCharges;
  depositionCharges.reserve(NDepositions);
  float totalDepositedCharge = 0;

  for (int i_dep = 0; i_dep < NDepositions; ++i_dep) {
    int charge = static_cast<int>(m_chargeSamplingHist->GetRandom(&randomGen));

    depositionCharges.push_back(charge);
    totalDepositedCharge += charge;
  }

  // rescale charges to ensure total charge is equal to simHit.charge()
  // will be important when drawing each depositions charge from the straggling distribution
  const float chargeScaler = simHit.charge() / static_cast<float>(totalDepositedCharge);
  for (int i_dep = 0; i_dep < NDepositions; ++i_dep) {
    float charge = depositionCharges[i_dep] * chargeScaler;

    // draw random position along the path for this deposition
    float t = randomGen.Rndm(); // uniform in (0,1)
    if (t < 0.f || t > 1.f)
      throw std::runtime_error("ChargeCollector_LUT::FillHit: Random position along path is out of bounds: t = " + std::to_string(t) + ". Must be in [0,1].");
    dd4hep::rec::Vector3D depositionPos = path.entry + t * path.travel;
    Index_voxel voxel;
    voxel.i = Trafo_local_pixIndex(depositionPos, m_digitizer.PixelPitch(), m_digitizer.PixelCount());
    voxel.j = Trafo_local_inpixIndex(depositionPos, m_LUT.GetBinCount(), m_digitizer.PixelPitch(), m_digitizer.ActiveVolumeDimensions());

    DistributeVoxelCharge(hitMap, voxel, charge, simHit);
  }

  m_digitizer.FillHistograms_fromChargeCollector_perSimHit(simHit.layer(), path.travel, path.lengthG4, simHit.truthPos(), trafoMatrix);
}

// void ChargeCollector_LUT::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix, TRandom3& randomGen) const {

//   Path path;
//   if (!ConstructPath(path, simHit, trafoMatrix, m_digitizer)) [[unlikely]]
//     return;

//   if (m_shiftTruthPos) {
//     MoveTruthPosition(simHit, path); // shifts the sim hit position to the depth in the sensor where most charge is collected, to get usesful residual plots.
//   }

//   const Index_inPix binCount = m_LUT.GetBinCount();
//   const std::array<float, 3> entry = {static_cast<float>(path.entry.x()), static_cast<float>(path.entry.y()), static_cast<float>(path.entry.z())};
//   const std::array<float, 3> travel = {static_cast<float>(path.travel.x()), static_cast<float>(path.travel.y()), static_cast<float>(path.travel.z())};

//   std::array<int, 3> g; // fine-grid bin index per axis, of the voxel the path currently is in
//   std::array<int, 3> step; // direction (+1/-1) the bin index moves along each axis
//   std::array<float, 3> tMax; // t at which the path crosses the next bin boundary on each axis
//   std::array<float, 3> tDelta; // t needed to cross one full bin on each axis
//   int iterationsLeft = 1; // upper bound on voxels crossed; guards against float quirks causing an endless loop

//   for (int ax = 0; ax < 3; ++ax) {
//     g[ax] = std::clamp(static_cast<int>(std::floor((entry[ax] - m_gridOrigin[ax]) / m_cellSize[ax])), 0, m_gridBinCount[ax] - 1);

//     if (travel[ax] != 0.f) {
//       step[ax] = (travel[ax] > 0.f) ? 1 : -1;
//       tDelta[ax] = m_cellSize[ax] / std::abs(travel[ax]);
//       const float nextBoundary = m_gridOrigin[ax] + (g[ax] + (step[ax] > 0 ? 1 : 0)) * m_cellSize[ax];
//       tMax[ax] = std::max(0.f, (nextBoundary - entry[ax]) / travel[ax]); // clamp to 0: the index clamping above can put the first boundary marginally behind the entry point
//       iterationsLeft += static_cast<int>(std::abs(travel[ax]) / m_cellSize[ax]) + 2;
//     }
//     else {
//       step[ax] = 0;
//       tDelta[ax] = std::numeric_limits<float>::infinity();
//       tMax[ax] = std::numeric_limits<float>::infinity();
//     }
//   }

//   Index_voxel vox;
//   vox.i = {g[0] / binCount[0], g[1] / binCount[1]};
//   vox.j = {g[0] % binCount[0], g[1] % binCount[1], g[2]};

//   float tPrev = 0.f;
//   while (iterationsLeft-- > 0) {
//     const int ax = (tMax[0] < tMax[1]) ? (tMax[0] < tMax[2] ? 0 : 2) : (tMax[1] < tMax[2] ? 1 : 2); // axis of the nearest boundary crossing. Corner ties are broken arbitrarily: the "wrong" voxel is visited with zero chord length, ie. zero charge
//     const float tNext = std::min(tMax[ax], 1.f);

//     if (tNext > tPrev) // skip zero-length chords (from corner ties or a path starting exactly on a boundary)
//       DistributeVoxelCharge(hitMap, vox, (tNext - tPrev) * simHit.charge(), simHit);

//     if (tMax[ax] >= 1.f)
//       break; // path ends inside the current voxel

//     tPrev = tMax[ax];
//     tMax[ax] += tDelta[ax];

//     g[ax] += step[ax];
//     if (g[ax] < 0 || g[ax] >= m_gridBinCount[ax]) [[unlikely]] {
//       if (1.f - tPrev > 1.e-4f)
//         m_digitizer.warning() << "ChargeCollector_LUT::FillHit: path left the voxel grid with a path fraction of " << 1.f - tPrev << " remaining. Dropping the corresponding charge." << endmsg;
//       break;
//     }

//     if (ax == 2) {
//       vox.j[2] += step[2]; // w has no pixel index; g range check above keeps j[2] valid
//     }
//     else {
//       vox.j[ax] += step[ax];
//       int& pixelIndex = (ax == 0) ? vox.i[0] : vox.i[1];
//       if (vox.j[ax] == binCount[ax]) {
//         vox.j[ax] = 0;
//         ++pixelIndex;
//       }
//       else if (vox.j[ax] < 0) {
//         vox.j[ax] = binCount[ax] - 1;
//         --pixelIndex;
//       }
//     }
//   } // voxel traversal

//   if (iterationsLeft < 0) [[unlikely]]
//     m_digitizer.warning() << "ChargeCollector_LUT::FillHit: voxel traversal did not terminate within the expected number of steps. Some charge may have been dropped." << endmsg;

//   m_digitizer.FillHistograms_fromChargeCollector_perSimHit(simHit.layer(), path.travel, path.lengthG4, simHit.truthPos(), trafoMatrix); // fill histograms once per sim hit, with info from the path (eg. travel vector, which contains info on the angle of incidence)
// }

void ChargeCollector_LUT::DistributeVoxelCharge(HitMap& hitMap, const Index_voxel& i_vox, const float charge, const SimHitWrapper& simHit) const {

  /* cache things, this is the hottest loop */
  const int lutSize = m_LUT.GetSize();
  const int i_u_origin = i_vox.i[0] - m_LUT.GetSizeHalf(); // pix index of leftmost pixel in LUT matrix
  const int i_v_origin = i_vox.i[1] - m_LUT.GetSizeHalf();
  const int pixelCount_u = static_cast<int>(m_digitizer.PixelCount().at(0));
  const int pixelCount_v = static_cast<int>(m_digitizer.PixelCount().at(1));

  const int col_min = std::max(0, -i_u_origin);
  const int col_max = std::min(lutSize, pixelCount_u - i_u_origin);

  const int row_min = std::max(0, -i_v_origin);
  const int row_max = std::min(lutSize, pixelCount_v - i_v_origin);

  for (int col = col_min; col < col_max; ++col) {
    const int i_u = i_u_origin + col; // convert from col in [0, matrixSize) to pixel offset in [-matrixSize_half, matrixSize_half]. Note size_half = (size-1)/2

    for (int row = row_min; row < row_max; ++row) {
      const int i_v = i_v_origin + row;

      const float chargeToAdd = m_LUT.GetWeight(i_vox.j, col, row) * charge;
      hitMap.FillCharge({i_u, i_v}, chargeToAdd, simHit);
    }
  }
}

void ChargeCollector_LUT::MoveTruthPosition(const SimHitWrapper& simHit, const Path& path) const {
  if (m_chargeCollectionDepthCenter == 0.f)
    return;

  /* shift the sim hit position along the path to the depth that is closest to the target w (ie. target depth) */

  float t = (m_chargeCollectionDepthCenter - path.entry.z()) / path.travel.z(); // how far along the path do we need to go to get to the target depth?
  // t in [0,1] means it's within the path, otherwise it's outside of the path and we will shift to the closest end (entry or exit)

  if (t < 0.f)
    t = 0.f;
  else if (t > 1.f)
    t = 1.f;

  simHit.SetTruthPos(path.entry + t * path.travel);
}

/* -- Single pixel approach -- */

ChargeCollector_SinglePixel::ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_SinglePixel constructed." << endmsg;
}

void ChargeCollector_SinglePixel::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix, TRandom3& randomGen) const {
  (void) trafoMatrix; // Not used in this implementation of ChargeCollector, but we need to keep it as argument to conform to the interface. Silences the unused parameter warning.
  (void) randomGen;

  const std::array<int, 2> i_uv = Trafo_local_pixIndex(simHit.truthPos(), m_digitizer.PixelPitch(), m_digitizer.PixelCount());

  if (!(i_uv[0] == -1 || i_uv[1] == -1 || i_uv[0] >= static_cast<int>(m_digitizer.PixelCount()[0]) || i_uv[1] >= static_cast<int>(m_digitizer.PixelCount()[1])))
    hitMap.FillCharge(i_uv, simHit.charge(), simHit);
}

/* -- Debug approach -- */

ChargeCollector_Debug::ChargeCollector_Debug(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_Debug constructed." << endmsg;
}

void ChargeCollector_Debug::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix, TRandom3& randomGen) const {
  (void) trafoMatrix; // Not used in this implementation of ChargeCollector, but we need to keep it as argument to conform to the interface. Silences the unused parameter warning.
  (void) randomGen;

  const dd4hep::rec::Vector3D pos_local = simHit.truthPos();
  const float charge = simHit.charge();
  const std::array<int, 2> i_uv = Trafo_local_pixIndex(pos_local, m_digitizer.PixelPitch(), m_digitizer.PixelCount());

  m_digitizer.verbose() << "     - SimHit at local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ")" << endmsg;
  m_digitizer.verbose() << "       - and pixel indices      (" << i_uv[0] << ", " << i_uv[1] << ")" << endmsg;
  if (i_uv[0] == -1 || i_uv[1] == -1 || i_uv[0] >= static_cast<int>(m_digitizer.PixelCount()[0]) || i_uv[1] >= static_cast<int>(m_digitizer.PixelCount()[1])) {
    m_digitizer.warning() << "simHit local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ") is out of sensor bounds U: [" << -m_digitizer.ActiveVolumeDimensions().at(0)/2 << ", " << m_digitizer.ActiveVolumeDimensions().at(0)/2 << "], V: [" << -m_digitizer.ActiveVolumeDimensions().at(1)/2 << ", " << m_digitizer.ActiveVolumeDimensions().at(1)/2 << "]. This simHit will be skipped." << endmsg;
  }
  else {
    m_digitizer.verbose() << "       - Filling charge " << simHit.charge() << " e." << endmsg;

    hitMap.FillCharge(i_uv, 0.5*charge, simHit);
    if (i_uv[0] + 1 < static_cast<int>(m_digitizer.PixelCount()[0]))
      hitMap.FillCharge({i_uv[0] + 1, i_uv[1]}, 0.3*charge, simHit);
    if (i_uv[1] + 1 < static_cast<int>(m_digitizer.PixelCount()[1]))
      hitMap.FillCharge({i_uv[0], i_uv[1] + 1}, 0.1*charge, simHit);
    if (i_uv[1] + 2 < static_cast<int>(m_digitizer.PixelCount()[1]))
      hitMap.FillCharge({i_uv[0], i_uv[1] + 2}, 0.1*charge, simHit);

    m_digitizer.verbose() << "       - Total charge collected in hitMap: " << hitMap.GetTotalCharge() << " e." << endmsg;
  }
}

// /* -- Drift approach -- */

} // n amespace VTXdigi_tools
