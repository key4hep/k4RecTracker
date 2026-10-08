// VTXdigi_Modular/src/ChargeCollector_impl.cpp

#include "../src/ChargeCollector_impl.h"
#include "VTXdigi_tools.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>

namespace VTXdigi_tools {
using ::VTXdigi_Modular; // "unqualified name introduction from global namespace" (just so I remember what to call this in C++ speak)

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

Path::Path(const SimHitWrapper& simHit, const TGeoHMatrix& trafoMatrix, const VTXdigi_Modular&  digitizer) {
  const bool printDebug = digitizer.msgLevel() <= MSG::DEBUG;
  simPos = simHit.truthPos();

  const double momentumEps = 1e-10;
  double momentum_global[3] = {
    static_cast<double>(simHit.hitPtr()->getMomentum().x),
    static_cast<double>(simHit.hitPtr()->getMomentum().y),
    static_cast<double>(simHit.hitPtr()->getMomentum().z)
  };
  if (std::abs(momentum_global[0]) < momentumEps && std::abs(momentum_global[1]) < momentumEps && std::abs(momentum_global[2]) < momentumEps) {
    digitizer.warning() << "SimHit momentum is zero (ie. smaller than " << momentumEps << "). Skipping simHit." << endmsg;
    isValid = false;
    return;
  }

  double momentum_local[3];
  trafoMatrix.MasterToLocalVect(momentum_global, momentum_local);
  const double momentum_local_normSq = momentum_local[0]*momentum_local[0] + momentum_local[1]*momentum_local[1] + momentum_local[2]*momentum_local[2];

  /* Compute linear approximation of the path */
  lengthG4 = simHit.hitPtr()->getPathLength();

  /* Step 1 - compute entry point & travel vector*/
  if (momentum_local[2]*momentum_local[2] / momentum_local_normSq < kMinPathCosTheta*kMinPathCosTheta) {
    digitizer.debug() << "SimHit momentum is (almost) parallel to the sensor surface (" << momentum_local[0] << ", " << momentum_local[1] << ", " << momentum_local[2] << ")." << endmsg;

    travel = dd4hep::rec::Vector3D(momentum_local[0], momentum_local[1], 0.);
    travel = 2. * lengthG4 * travel.unit(); // double the length to make sure the clipping in step 2 does not cut off too much (eg. when the pos is exactly at the sensor edge)
    entry = simPos - 0.5 * travel;
  }
  else{
    const double scaleFactor_travel = digitizer.ActiveVolumeDimensions().at(2) / std::abs(momentum_local[2]);
    travel = scaleFactor_travel * dd4hep::rec::Vector3D(momentum_local[0], momentum_local[1], momentum_local[2]);

    double shiftDist_w;
    if (travel.z() >= 0.0) {
      shiftDist_w = simPos.z() + 0.5 * digitizer.ActiveVolumeDimensions().at(2);
    }
    else {
      shiftDist_w = simPos.z() - 0.5 * digitizer.ActiveVolumeDimensions().at(2);
    }
    const double scaleFactor_entry = shiftDist_w / travel.z();
    entry = simPos - scaleFactor_entry * travel;
  }

  /* Step 2 - clip path to sensor edges (in u/v) */
  std::array<double, 2> t = {0., 1.}; // parametrize path as entry + t*travel; t in [0,1]
  t = ComputePathClippingFactors(t, entry.x(), travel.x(), digitizer.ActiveVolumeDimensions().at(0));
  t = ComputePathClippingFactors(t, entry.y(), travel.y(), digitizer.ActiveVolumeDimensions().at(1));
  if (t[0] != 0.0 || t[1] != 1.0) {
    if (0.0 <= t[0] && t[0] < t[1] && t[1] <= 1.0) {
      /* valid clipping */
      if (printDebug) digitizer.debug() << "       - Clipping SimHitPath with t [" << t[0] << ", " << t[1] << "]. PathLength changed to " << static_cast<int>((t[1] - t[0]) * travel.r()*1000) << " um from " << static_cast<int>(travel.r()*1000) << " um" << endmsg;

      entry = entry + t[0] * travel;
      travel = (t[1] - t[0]) * travel;
    }
    else {
      /* invalid clipping, shouldn't happen */
      digitizer.warning() << "VTXdigi_tools::ConstructPath() - invalid clipping factors t = [" << t[0] << ", " << t[1] << "]. Path might lie completely outside the sensor." << endmsg;
      digitizer.debug() << " -> entry (" << entry.x() << ", " << entry.y() << ", " << entry.z() << ") mm, exit (" << entry.x() + travel.x() << ", " << entry.y() + travel.y() << ", " << entry.z() + travel.z() << ") mm, sensor dim. (+-" << digitizer.ActiveVolumeDimensions().at(0)/2 << ", +-" << digitizer.ActiveVolumeDimensions().at(1)/2 << ") mm" << endmsg;
      digitizer.debug() << " -> Path length " << static_cast<int>(travel.r()*1000) << " um, in G4 " << static_cast<int>(simHit.hitPtr()->getPathLength()*1000) << " um" << endmsg;
      isValid = false;
      return;
    }
  }

  /* Step 3 -check that path is not much longer than the length it had in Geant4. Order of steps 2 and 3 is important! */
  if (travel.r() > kPathLengthTolerance * lengthG4) {
    if (printDebug) digitizer.debug() << "       - Shortening path length from " << static_cast<int>(travel.r()*1000) << " um to " << static_cast<int>(lengthG4*1000) << " um (the respective path length in Geant4)." << endmsg;

    /* make sure the path stays centred around the simTrackerHit position */
    const double t_simPos = ( (simPos - entry).dot(travel) ) / (travel.r() * travel.r());

    const double t_length_halved = 0.5 * lengthG4 / travel.r(); // length of the new path in terms of t [0,1] on old path, halved
    const double t_center = std::max(t_length_halved, std::min(t_simPos, 1. - t_length_halved)); // center of new path clamped to [t_length_half, 1 - t_length_half] while not exceeding [0,1]

    const double t_min = t_center - t_length_halved;
    const double t_max = t_center + t_length_halved;

    entry = entry + t_min * travel;
    travel = (t_max - t_min) * travel;
  }

  if (printDebug) digitizer.debug() << "       - Constructed path, length " << travel.r()*1000 << " um (G4-length " << lengthG4*1000 << " um), entry (" << entry.x() << ", " << entry.y() << ", " << entry.z() << ") mm, exit (" << entry.x() + travel.x() << ", " << entry.y() + travel.y() << ", " << entry.z() + travel.z() << ") mm, " << endmsg;
  isValid = true;
}

std::vector<std::pair<float, dd4hep::rec::Vector3D>> Path::SampleDepositions(const float hitCharge, TRandom3& randomGen, const float meanDepositionsPerUm, const TH1D& chargeSamplingHist) const {
  // draw number of depositions from a sub-poissonian distribution
  float width = 0.8; // width of the sub-poissonian (0.8 resembles what Allpix Squared does very well)
  float mean = travel.r() * 1000.f * meanDepositionsPerUm ; // convert from mm to um
  int NDepositions = std::max(1, static_cast<int>(randomGen.Binomial(std::round(mean / width), width)));

  std::vector<std::pair<float, dd4hep::rec::Vector3D>> depositions;
  depositions.reserve(NDepositions);
  float totalDepositedCharge = 0;

  for (int i_dep = 0; i_dep < NDepositions; ++i_dep) {
    int charge = static_cast<int>(chargeSamplingHist.GetRandom(&randomGen));

    const double t = randomGen.Rndm(); // uniform in (0,1)
    const dd4hep::rec::Vector3D pos = entry + t * travel;

    depositions.push_back(std::make_pair(charge, pos));
    totalDepositedCharge += charge;
  }

  // rescale charges to ensure total charge is equal to simHit.charge()
  // will be important when drawing each depositions charge from the straggling distribution
  const float chargeScaler = hitCharge / static_cast<float>(totalDepositedCharge);
  for (auto& dep : depositions) {
    dep.first *= chargeScaler;
  }
  return depositions;
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
  const bool printDebug = digitizer.msgLevel() <= MSG::DEBUG;

  if (printDebug) digitizer.debug() << " - Constructing LUT from file \"" << lutFileName << "\"." << endmsg;

  /* parse the LUT in Allpix Squared format
   * See https://indico.cern.ch/event/1489052/contributions/6475539/attachments/3063712/5418424/Allpix_workshop_Lemoine.pdf (slide 10) for more info on fields in the LUT file */

  if (lutFileName.empty())
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): LUT file name is empty. A LUT file must be given to load the lookup table.");

  const int headerLines = 5;

  if (printDebug) digitizer.debug() << "   - Opening LUT file \"" << lutFileName << "\"." << endmsg;
  std::ifstream lutFile(lutFileName);
  if (!lutFile.is_open())
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Could not open LUT file \"" + lutFileName + "\".");

  std::string line;
  int lineCount = 0;

  /* Parse pixel-pitch, thickness, in-pixel bin count from header (all in 5th line) */
  for (; lineCount < 5; ++lineCount)
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
    m_voxelCount.at(j) = std::stoi(headerLineEntries.at(7+j));
  }
  if (printDebug) digitizer.debug() << "   - found in-pixel bin count of (" << m_voxelCount.at(0) << ", " << m_voxelCount.at(1) << ", " << m_voxelCount.at(2) << ") from LUT file header." << endmsg;

  /* -> compare the values we just parsed to the values retrieved from the detector geometry */
  const double eps = 1e-12; // reasonable for number O(0.01) (like sensor thickness in mm) with double precision

  const double sensorThickness = std::stod(headerLineEntries.at(0)) / 1000.0; // convert from um to mm
  if (std::abs(sensorThickness - digitizer.ActiveVolumeDimensions().at(2)) > eps) {
    if (!digitizer.LUT_ignorePitch()) {
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Sensor thickness mismatch between LUT file and detector geometry: LUT file specifies " + std::to_string(sensorThickness) + " mm, but geometry has " + std::to_string(digitizer.ActiveVolumeDimensions().at(2)) + " mm active volume thickness.");
    }
    else {
      digitizer.warning() << "Sensor thickness mismatch between LUT file and detector geometry. LUT file: " << sensorThickness << "mm, geometry (active volume thickness): " << digitizer.ActiveVolumeDimensions().at(2) << "mm. Ignored because LookupTableIgnorePitch is set to true." << endmsg;
    }
  }

  const std::array<double, 2> pitch = {std::stod(headerLineEntries.at(1)) / 1000.0, std::stod(headerLineEntries.at(2)) / 1000.0};
  if (std::abs(pitch[0] - digitizer.PixelPitch().at(0)) > eps || std::abs(pitch[1] - digitizer.PixelPitch().at(1)) > eps) {
    if (!digitizer.LUT_ignorePitch())
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Pixel pitch mismatch between LUT file and detector geometry: LUT file specifies (" + std::to_string(pitch[0]) + ", " + std::to_string(pitch[1]) + ") mm, but geometry has (" + std::to_string(digitizer.PixelPitch().at(0)) + ", " + std::to_string(digitizer.PixelPitch().at(1)) + ") mm.");
    else
      digitizer.warning() << "Pixel pitch mismatch between LUT file and detector geometry. LUT file: " << pitch[0] << "mm, geometry: " << digitizer.PixelPitch().at(0) << "mm. Ignored because LookupTableIgnorePitch is set to true." << endmsg;
  }

  if (printDebug) digitizer.debug() << "   - Found matching pixel pitch and sensor thickness in LUT file." << endmsg;

  /* Get matrix size (5x5, 7x7, ...) from the length of the first line after the header */
  bool foundDataLine = false;
  while (std::getline(lutFile, line)) {
    if (line.empty() || line[0] == '#') {
      if (printDebug) digitizer.debug() << "VTXdigi_tools::LookupTable::LookupTable(): Empty or comment line found in LUT file at line " << lineCount+1 << ". Ignoring" << endmsg;
      continue;
    }

    // count whitespace-separated tokens, robust against repeated spaces, tabs and trailing whitespace / CR
    std::istringstream lineStream(line);
    std::string token;
    int tokenCount = 0;
    while (lineStream >> token)
      ++tokenCount;

    const int entryCount = tokenCount - 3; // first 3 entries are bin indices
    m_matrixSize = static_cast<int>(std::lround(std::sqrt(std::max(entryCount, 0))));
    if (m_matrixSize * m_matrixSize != entryCount)
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): First data line in LUT file has " + std::to_string(entryCount) + " matrix entries (after 3 bin indices), which is not a perfect square. File: " + digitizer.LutFileName());

    foundDataLine = true;
    break;
  }
  if (!foundDataLine)
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Could not find a data line after header in LUT file: " + digitizer.LutFileName());

  if (m_matrixSize < 3 || m_matrixSize % 2 == 0)
    throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Matrix size must be an odd integer >= 3, but is " + std::to_string(m_matrixSize) + ".");
  m_matrixSize_half = (m_matrixSize - 1) / 2;
  if (printDebug) digitizer.debug() << "   - Inferred matrix size of " << m_matrixSize << " from first line." << endmsg;

  /* Set up the matrix vector */
  m_matrices.resize(m_voxelCount.at(0) * m_voxelCount.at(1) * m_voxelCount.at(2) * m_matrixSize * m_matrixSize, 0.f);

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

  if (printDebug) digitizer.debug() << "   - Parsing LUT file, filling into lookup table." << endmsg;
  std::vector<float> matricesEntrySum_perWBin;
  matricesEntrySum_perWBin.resize(m_voxelCount.at(2), 0.f);
  float matricesEntrySum = 0.f;

  lineCount = headerLines + 1;
  while (std::getline(lutFile, line)) {
    if (line.empty() || line[0] == '#') {
      if (printDebug)  digitizer.debug() << "VTXdigi_tools::LookupTable::LookupTable(): Empty or comment line found in LUT file at line " << lineCount+1 << ". Ignoring" << endmsg;
      continue;
    }

    std::istringstream stringStream(line);
    std::vector<std::string> lineEntries;
    std::string entryString;

    /* read the line & do sanity checks */
    while (std::getline(stringStream, entryString, ' ')) {
      if (!entryString.empty())
        lineEntries.push_back(entryString);
    }

    if (static_cast<int>(lineEntries.size()) != 3 + m_matrixSize*m_matrixSize)
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Invalid number of entries in LUT file at line " + std::to_string(lineCount+1) + ": found " + std::to_string(lineEntries.size()) + " entries, but expected " + std::to_string(3 + m_matrixSize*m_matrixSize) + " (3 for bin indices, " + std::to_string(m_matrixSize*m_matrixSize) + " for matrix values).");

    /* First 3 entries are in-pixel binning indices */
    VoxelIndex voxI({std::stoi(lineEntries[0])-1, std::stoi(lineEntries[1])-1, std::stoi(lineEntries[2])-1});// Allpix2 input is 1-indexed. Insane, I know.

    if (voxI.at(0) < 0 || voxI.at(0) >= m_voxelCount[0] ||
        voxI.at(1) < 0 || voxI.at(1) >= m_voxelCount[1] ||
        voxI.at(2) < 0 || voxI.at(2) >= m_voxelCount[2]) {
      throw std::runtime_error("Invalid in-pixel bin indices in LUT file at line " + std::to_string(lineCount+1) + ": got (" + std::to_string(voxI.at(0)) + ", " + std::to_string(voxI.at(1)) + ", " + std::to_string(voxI.at(2)) + "), but expected ranges are [0, " + std::to_string(m_voxelCount[0]-1) + "], [0, " + std::to_string(m_voxelCount[1]-1) + "], [0, " + std::to_string(m_voxelCount[2]-1) + "].");
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
    // digitizer.verbose() << "   - Parsed matrix for in-pixel bin (" << voxI.at(0) << ", " << voxI.at(1) << ", " << voxI.at(2) << "), entry sum " << std::to_string(matrixEntrySum) << ", setting it now..." << endmsg;
    if (std::isnan(matrixEntrySum))
      throw std::runtime_error("VTXdigi_tools::LookupTable::LookupTable(): Charge sharing matrix for in-pixel bin (" + std::to_string(voxI.at(0)) + "," + std::to_string(voxI.at(1)) + "," + std::to_string(voxI.at(2)) + ") contains NaN values (sum of entries is NaN).");

    matricesEntrySum += matrixEntrySum;
    matricesEntrySum_perWBin[voxI.at(2)] += matrixEntrySum;
    SetMatrix(voxI, matrixEntries);

    lineCount++;
  } // loop over lines containing a matrix each

  const int matrixCount = lineCount - (headerLines+1);

  if (matrixCount != m_voxelCount[0] * m_voxelCount[1] * m_voxelCount[2])
    throw std::runtime_error("Invalid number of matrices loaded from file: expected " + std::to_string(m_voxelCount[0] * m_voxelCount[1] * m_voxelCount[2]) + " matrices (inferred from bin count in header) but found " + std::to_string(matrixCount) + " lines.");

  // From which w-level are charges collected?
  double collectedFromW = 0.0;
  for (int j_w = 0; j_w < m_voxelCount.at(2); j_w++) {
    const double binHeight = sensorThickness / m_voxelCount.at(2);
    const double w_bin_center = j_w*binHeight + 0.5*binHeight - 0.5*sensorThickness;
    collectedFromW += matricesEntrySum_perWBin.at(j_w) * w_bin_center;
  }
  collectedFromW /= matricesEntrySum;

  if (std::abs(collectedFromW) > 0.5 * sensorThickness) throw std::runtime_error("Invalid charge collection depth inferred from LUT file: " + std::to_string(collectedFromW) + " mm. This is outside the sensor volume, which extends from " + std::to_string(-0.5*sensorThickness) + " mm to " + std::to_string(0.5*sensorThickness) + " mm.");
  if (std::abs(collectedFromW) >= 1.e-5) {
    m_chargeCollectionDepthCenter = static_cast<double>(collectedFromW); // the member is initalised to 0.0
  }

  digitizer.info() << " - Loaded lookup table from file. Matrices parsed: " << (matrixCount) << ". Charge collected from sensitive volume: " << matricesEntrySum/static_cast<double>(matrixCount)*100 << " percent (rest is lost, eg recombination). The center of the collection region is at " << m_chargeCollectionDepthCenter << endmsg;
}

void LookupTable::SetMatrix(const VoxelIndex& voxI, const std::vector<float>& weights) {
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
    throw std::runtime_error("VTXdigi_tools::LookupTable::SetMatrix: Charge sharing matrix for in-pixel bin (" + std::to_string(voxI.at(0)) + "," + std::to_string(voxI.at(1)) + "," + std::to_string(voxI.at(2)) + ") contains NaN values.");
  if (sum < 0 || sum > 1.f + 1.e-5f)
    throw std::runtime_error("VTXdigi_tools::LookupTable::SetMatrix: Charge sharing matrix for in-pixel bin (" + std::to_string(voxI.at(0)) + "," + std::to_string(voxI.at(1)) + "," + std::to_string(voxI.at(2)) + ") has a weight sum of " + std::to_string(sum) + ", but needs to lie in [0,1].");

  for (int row = 0; row < m_matrixSize; ++row) {
    for (int col = 0; col < m_matrixSize; ++col) {
      /* weights are given in row-major order, starting at top left.
        * We store charge sharing matrices in col-major order, starting at bottom left (lowest bin index) */
      m_matrices.at(FindIndex(voxI, col, row)) = weights.at((m_matrixSize-1-row)*m_matrixSize + col);
    }
  }
}

void LookupTable::SetAllMatrices(const std::vector<float>& weights) {
  for (int j_u = 0; j_u < m_voxelCount.at(0); ++j_u) {
    for (int j_v = 0; j_v < m_voxelCount.at(1); ++j_v) {
      for (int j_w = 0; j_w < m_voxelCount.at(2); ++j_w) {
        SetMatrix({j_u, j_v, j_w}, weights);
      }
    }
  }
}

int LookupTable::FindIndex (const VoxelIndex& voxI, const int col, const int row) const {
  #ifndef NDEBUG
    if (voxI[0] < 0 || voxI[0] >= m_voxelCount[0]
      || voxI[1] < 0 || voxI[1] >= m_voxelCount[1]
      || voxI[2] < 0 || voxI[2] >= m_voxelCount[2] ) {
      throw std::runtime_error("VTXdigi_tools::LookupTable::FindIndex: in-pix bin out of range");
    }
    if (col < 0 || col >= m_matrixSize || row < 0 || row >= m_matrixSize) {
      throw std::runtime_error("VTXdigi_tools::LookupTable::FindIndex: col or row out of range");
    }
  #endif

  int index_matrix = voxI[0] + m_voxelCount[0] * (voxI[1] + m_voxelCount[1] * voxI[2]);
  int index_element = col * m_matrixSize + row;
  return index_matrix * m_matrixSize * m_matrixSize + index_element;
}

ChargeCollector_LUT::ChargeCollector_LUT(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer),
  m_LUT(digitizer.LutFileName(), digitizer),
  m_shiftTruthPos(digitizer.LUT_shiftTruthPos()) {

  /* LUT is constructed in place (from file) */
  m_chargeCollectionDepthCenter = m_LUT.GetChargeCollectionDepthCenter();

  /* Load charge deposition sampling parameters */
  m_meanDepositionsPerUm = digitizer.MeanDepositionsPerUm();

  // shipped with k4RecTracker, VTXDIGI_MODULAR_DATADIR is set at compile time in CMakeLists.txt
  const std::string chargeDepFileName = VTXDIGI_MODULAR_DATADIR "/chargeDepositionDistribution.root";
  m_digitizer.info() << " - Loading deposition charge histogram from " << chargeDepFileName << endmsg;

  std::unique_ptr<TFile> file(TFile::Open(chargeDepFileName.c_str(), "READ"));
  if (!file || file->IsZombie())
    throw std::runtime_error("Could not open deposition charge histogram file " + chargeDepFileName + ", cannot continue.");

  TH1D* hist_chargeDep = file->Get<TH1D>("deposition_charge");
  if (!hist_chargeDep)
    throw std::runtime_error("Could not find histogram \"deposition_charge\" in file "+ chargeDepFileName + ", cannot continue.");

  // move ownership from TFile to this class
  hist_chargeDep->SetDirectory(nullptr);
  hist_chargeDep->ComputeIntegral(); // so that we don't compute the integral later down the line (which would be slow & I am not sure about thread safety)
  m_chargeSamplingHist.reset(hist_chargeDep);

  m_digitizer.info() << " - ChargeCollector_LUT constructed successfully." << endmsg;
}

void ChargeCollector_LUT::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix, TRandom3& randomGen) const {

  Path path(simHit, trafoMatrix, m_digitizer);
  if (!path.isValid) [[unlikely]]
    return;

  if (m_shiftTruthPos) {
    MoveTruthPosition(simHit, path); // shifts the sim hit position to the depth in the sensor where most charge is collected, to get usesful residual plots.
  }

  const auto pixelPitch = m_digitizer.PixelPitch();
  const auto pixelCount = m_digitizer.PixelCount();
  const auto thickness = m_digitizer.ActiveVolumeDimensions().at(2);
  const auto lutVoxelCount = m_LUT.GetVoxelCount();
  for (const auto& [charge, pos] : path.SampleDepositions(simHit.charge(), randomGen, m_meanDepositionsPerUm, *m_chargeSamplingHist)) {
    const auto [pixI, voxI] = VTXdigi_tools::Trafo_local_pixIVoxI(pos, pixelPitch, pixelCount, thickness, lutVoxelCount);

    DistributeVoxelCharge(hitMap, pixI, voxI, charge, simHit);
  }

  m_digitizer.FillHistograms_fromChargeCollector_perSimHit(simHit.layer(), path.travel, path.lengthG4, simHit.truthPos(), trafoMatrix);
}

void ChargeCollector_LUT::DistributeVoxelCharge(HitMap& hitMap, const PixelIndex& pixI, const VoxelIndex& voxI, const float charge, const SimHitWrapper& simHit) const {

  /* cache things, this is the hottest loop */
  const int lutSize = m_LUT.GetSize();
  const int i_u_origin = pixI[0] - m_LUT.GetSizeHalf(); // pix index of leftmost pixel in LUT matrix
  const int i_v_origin = pixI[1] - m_LUT.GetSizeHalf();
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

      const float chargeToAdd = m_LUT.GetWeight(voxI, col, row) * charge;
      hitMap.FillCharge({i_u, i_v}, chargeToAdd, simHit);
    }
  }
}

void ChargeCollector_LUT::MoveTruthPosition(const SimHitWrapper& simHit, const Path& path) const {
  if (m_chargeCollectionDepthCenter == 0.0)
    return;

  /* shift the sim hit position along the path to the depth that is closest to the target w (ie. target depth) */

  double t = (m_chargeCollectionDepthCenter - path.entry.z()) / path.travel.z(); // how far along the path do we need to go to get to the target depth?
  // t in [0,1] means it's within the path, otherwise it's outside of the path and we will shift to the closest end (entry or exit)

  if (t < 0.0)
    t = 0.0;
  else if (t > 1.0)
    t = 1.0;

  simHit.SetTruthPos(path.entry + t * path.travel);
}

/* -- Single pixel approach -- */

ChargeCollector_SinglePixel::ChargeCollector_SinglePixel(const VTXdigi_Modular& digitizer) : IChargeCollector(digitizer) {
  m_digitizer.debug() << "ChargeCollector_SinglePixel constructed." << endmsg;
}

void ChargeCollector_SinglePixel::FillHit(const SimHitWrapper& simHit, HitMap& hitMap, const TGeoHMatrix& trafoMatrix, TRandom3& randomGen) const {
  (void) trafoMatrix; // Not used in this implementation of ChargeCollector, but we need to keep it as argument to conform to the interface. Silences the unused parameter warning.
  (void) randomGen;

  const PixelIndex pixI = Trafo_local_pixI(simHit.truthPos(), m_digitizer.PixelPitch(), m_digitizer.PixelCount());
  hitMap.FillCharge(pixI, simHit.charge(), simHit);
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
  const PixelIndex pixI = Trafo_local_pixI(pos_local, m_digitizer.PixelPitch(), m_digitizer.PixelCount());





  hitMap.FillCharge(pixI, 0.5*charge, simHit);
  if (pixI[0] + 1 < static_cast<int>(m_digitizer.PixelCount()[0]))
    hitMap.FillCharge({pixI[0] + 1, pixI[1]}, 0.3*charge, simHit);
  if (pixI[1] + 1 < static_cast<int>(m_digitizer.PixelCount()[1]))
    hitMap.FillCharge({pixI[0], pixI[1] + 1}, 0.1*charge, simHit);
  if (pixI[1] + 2 < static_cast<int>(m_digitizer.PixelCount()[1]))
    hitMap.FillCharge({pixI[0], pixI[1] + 2}, 0.1*charge, simHit);

  if (m_digitizer.msgLevel() <= MSG::VERBOSE) {
    m_digitizer.verbose() << "     - Filling pixels for SimHit at local position (" << pos_local.x() << ", " << pos_local.y() << ", " << pos_local.z() << ")" << endmsg;
    m_digitizer.verbose() << "       - and pixel indices                         (" << pixI[0] << ", " << pixI[1] << ")" << endmsg;
    m_digitizer.verbose() << "       - Charge " << simHit.charge() << " e. Total charge collected in hitMap: " << hitMap.GetTotalCharge() << " e." << endmsg;
  }
}

} // n amespace VTXdigi_tools
