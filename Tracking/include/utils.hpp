#ifndef UTILS_HPP
#define UTILS_HPP

//=== Standard Library ===
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

//=== DD4hep / DDRec ===
#include "DD4hep/DetElement.h"
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include <DDRec/DetectorData.h>

//=== edm4hep ===
#include "edm4hep/TrackState.h"
#include "extension/MutableTrack.h"
#include "extension/TrackCollection.h"

#include <Gaudi/Algorithm.h>

//=== Others ===
#include <Objects/Helix.h>
#include <TMatrixDSym.h>
#include <TVector3.h>
#include <TVectorD.h>
#include <marlinutil/HelixClass_double.h>
#include <torch/torch.h>

/////////////////////////////////
/// Track propagation to Ecal ///
/////////////////////////////////

/**
 * @brief Retrieves a unique LayeredCalorimeterData extension from the detector.
 *
 * This function queries the global DD4hep detector instance for DetElements
 * matching the given inclusion and exclusion flags. It expects exactly one
 * DetElement to be selected. If none or multiple DetElements are found, it
 * throws a std::runtime_error with an appropriate error message.
 *
 * Once the unique DetElement is found, its extension of type
 * dd4hep::rec::LayeredCalorimeterData is returned.
 *
 * @param includeFlag Bitmask flag(s) specifying which detector types to include.
 * @param excludeFlag Bitmask flag(s) specifying which detector types to exclude.
 *
 * @return Pointer to the dd4hep::rec::LayeredCalorimeterData extension of the selected DetElement.
 *
 * @throws std::runtime_error If the global detector instance is not initialized.
 * @throws std::runtime_error If the selection of DetElements is not unique or empty.
 */
dd4hep::rec::LayeredCalorimeterData* getExtension(unsigned int includeFlag, unsigned int excludeFlag);

/**
 * @brief Computes the track state extrapolated at the calorimeter surface.
 *
 * This function takes the projected position of the track at the calorimeter
 * (ecalProjection), the helix representation of the track at its last hit,
 * and the magnetic field along the z-axis (Bz). It extrapolates the track
 * momentum to the calorimeter position and creates a new helix at this point.
 *
 * The output is an edm4hep::TrackState filled with the helix parameters at the
 * calorimeter, including position, momentum direction, curvature, and reference point.
 *
 * @param ecalProjection CartesianVector representing the track position at the calorimeter.
 * @param helixAtLastHit HelixClass_double instance describing the track at the last hit point.
 * @param Bz Magnetic field value along the z-axis in the detector region.
 *
 * @return edm4hep::TrackState The extrapolated track state at the calorimeter surface.
 */
edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection,
                                                  const HelixClass_double& helixAtLastHit, double m_Bz);

/**
 * @brief Extrapolates the track parameters to the calorimeter surface and updates the track with this information.
 *
 * This function takes a mutable track object and magnetic field parameters to
 * extrapolate the track helix from the last measured hit position to the
 * calorimeter surface (either barrel or endcap). It calculates the track
 * momentum at the last hit, constructs a helix representation, and then
 * projects the track onto the calorimeter surfaces to find the best
 * extrapolation point.
 *
 * The extrapolation considers both the barrel and endcap calorimeter
 * geometries and selects the projection with the earliest arrival time.
 * The corresponding TrackState at the calorimeter is then attached to the
 * provided track object.
 *
 * @param edm4hep_track MutableTrack object to be updated with calorimeter extrapolation.
 * @param m_Bz Magnetic field component along the z-axis.
 * @param charge Electric charge of the particle track.
 * @param a Conversion constant relating curvature to momentum (e.g., 0.299792458 GeV/(T*m)).
 * @param m_eCalBarrelInnerR Inner radius of the barrel calorimeter.
 * @param m_eCalBarrelMaxZ Maximum half-length (z) of the barrel calorimeter.
 * @param m_eCalEndCapInnerR Inner radius of the endcap calorimeter.
 * @param m_eCalEndCapOuterR Outer radius of the endcap calorimeter.
 * @param m_eCalEndCapInnerZ z-position of the inner surface of the endcap calorimeter.
 */
void FillTrackWithCalorimeterExtrapolation(extension::MutableTrack& edm4hep_track, double m_Bz, int charge, double a,
                                           double m_eCalBarrelInnerR, double m_eCalBarrelMaxZ,
                                           double m_eCalEndCapInnerR, double m_eCalEndCapOuterR,
                                           double m_eCalEndCapInnerZ);

////////////////////////////
/// Clustering Utilities ///
////////////////////////////

/**
 * @brief Finds indices of "conditional points" based on a threshold and selection mask.
 *
 * This function identifies points from the `betas` tensor that exceed the given
 * threshold `tbeta` and are also present in the `unassigned` set. It returns
 * the indices of these points sorted by the negative value of their beta, effectively
 * sorting them in ascending order of beta.
 *
 * @param betas Tensor of beta values for all points (shape: [size_b]).
 * @param unassigned Tensor containing indices of points that are currently unassigned.
 * @param tbeta Threshold value for beta to select conditional points.
 *
 * @return Tensor of indices of conditional points sorted by ascending beta value.
 */
torch::Tensor find_condpoints(const torch::Tensor& betas, const torch::Tensor& unassigned, float tbeta);

/**
 * @brief Performs clustering on a set of points using beta thresholding and distance criteria.
 *
 * Given an input vector representing points with associated beta values,
 * this function clusters points based on a beta threshold (`tbeta`) and a
 * maximum distance (`td`). It identifies conditional points where beta
 * exceeds `tbeta`, then iteratively assigns points within distance `td` of
 * each conditional point to the same cluster.
 *
 * @param output_vector Vector of floats containing point data in the format [x, y, z, beta].
 * @param num_rows Number of points (rows) in the input data.
 * @param tbeta Beta threshold for selecting conditional points.
 * @param td Distance threshold for clustering points around conditional points.
 *
 * @return A 1D torch::Tensor of type Long containing cluster assignments for each point.
 *         Points assigned to cluster 0 are considered unclustered.
 */
torch::Tensor get_clustering(const std::vector<float>& output_vector, int num_rows, float tbeta, float td);

///////////////////////////////////////
/// Miscellaneous Utility Functions ///
///////////////////////////////////////

/**
 * @brief Returns the hypothesized charge for a given PDG code.
 *
 * Maps common particle PDG codes to their charge:
 * - Electrons (11): -1
 * - Positrons (-11): +1
 * - Muons (13): -1
 * - Anti-muons (-13): +1
 * - Charged pions (211, -211): ±1
 * - Charged kaons (321, -321): ±1
 * - Protons (2212): +1
 * - Anti-protons (-2212): -1
 *
 * Returns 0 if the PDG code is not recognized.
 *
 * @param pdg Particle Data Group (PDG) code of the particle.
 * @return int Hypothesized electric charge of the particle.
 */
int getHypotesisCharge(int pdg);

/**
 * @brief Computes the covariance matrix of the track state parameters.
 *
 * The state parameters are converted from Cartesian coordinates (position, momentum, time)
 * to track parameters (omega, phi0, d0, z0, tanLambda, time), and their covariance is propagated
 * through the Jacobian matrix J.
 *
 * Units expected:
 * - positions in mm
 * - momentum in GeV
 * - time in ns
 *
 * @param stateTrack 6D state vector: (x, y, z, px, py, pz)
 * @param params 6D vector of track parameters: (omega, phi0, d0, z0, tanLambda, time)
 * @param referencePoint 3D point for impact parameter calculation - see
 * https://flc.desy.de/lcnotes/notes/localfsExplorer_read?currentPath=/afs/desy.de/group/flc/lcnotes/LC-DET-2006-004.pdf
 * @param timeError uncertainty on the time parameter (ns)
 * @param statecovMatrix covariance matrix of the original state vector (6x6)
 * @return TMatrixDSym covariance matrix of the track parameters + time (6x6)
 */
TMatrixDSym computeTrackStateCovMatrix(TVectorD stateTrack, TVectorD params, TVector3 referencePoint, double timeError,
                                       TMatrixDSym statecovMatrix);

/**
 * @brief Print the message in a block of stars ('*') at level DEBUG
 */
void printInStars(const Gaudi::Algorithm* thisAlg, const std::string& msg, const int lineWidth);

#endif // UTILS_HPP
