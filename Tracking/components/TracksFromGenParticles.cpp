#include "Gaudi/Property.h"

// edm4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"

// marlin
#include <marlinutil/HelixClass_double.h>

// pandora
#include <Objects/Helix.h>

// k4FWCore
#include "k4FWCore/Transformer.h"

// DD4HEP
#include <DDRec/DetectorData.h>
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Readout.h"

// C++
#include <string>

using SimTrackerHitColl = std::vector<const edm4hep::SimTrackerHitCollection*>;

/** @class TracksFromGenParticles
 *
 *  Gaudi transformer that builds an edm4hep::TrackCollection out of an edm4hep::MCParticleCollection.
 *  It just builds an helix out of the genParticle position, momentum, charge and z component of the (constant) magnetic field, retrieved from the detector.
 *  From this helix, different edm4hep::TrackStates (AtIP, AtFirstHit, AtLastHit) are defined.
 *  The first and last hits are defined as those with smallest and largest time in the input SimTrackerHit collections
 *  The algorithm also performs extrapolation to the EM calorimeter inner face. This is done using the positions of the barrel and endcap retrieved from the detector data extensions.
 *  This is meant to enable technical development needing edm4hep::Track and performance studies where having generator based tracks is a reasonable approximation.
 *  GM TODO:
 *  - we could replace the generator-level SimTrackerHit collections with the digitised TrackerHit3D collections and create associations between the "reconstructed"
 *    tracks and the digitised hits (how to associate hit<->track? hit->sim hit->particle? or geometric matching?)
 *  - in case of intersections of the extrapolation with both ECAL barrel and endcap inner faces, we could keep both (could be useful for reconstruction)
 *    instead of only the one with lower arrival time
 *  @author Brieuc Francois
 *  @author Archil Durglishvili
 *  @author Giovanni Marchiori
 */

struct TracksFromGenParticles final
  : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>(const edm4hep::MCParticleCollection&, const SimTrackerHitColl&)> {

  TracksFromGenParticles(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(
            name, svcLoc,
            {KeyValues("InputGenParticles", {"MCParticles"}),
             KeyValues("InputSimTrackerHits", {"SimTrackerHits"})},
            {KeyValues("OutputTracks", {"TracksFromGenParticles"}),
             KeyValues("OutputMCRecoTrackParticleAssociation", {"TracksFromGenParticlesAssociation"})}) {
  }

  double getFieldFromCompact() {
    dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
    const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
    double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
    mainDetector.field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep
    return magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
  }

  dd4hep::rec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag) {

    dd4hep::rec::LayeredCalorimeterData * theExtension = 0;

    dd4hep::Detector & mainDetector = dd4hep::Detector::getInstance();
    const std::vector< dd4hep::DetElement>& theDetectors = dd4hep::DetectorSelector(mainDetector).detectors(  includeFlag, excludeFlag );

    debug() << " getExtension :  includeFlag: " << dd4hep::DetType( includeFlag ) << " excludeFlag: " << dd4hep::DetType( excludeFlag )
            << "  found : " << theDetectors.size() << "  - first det: " << theDetectors.at(0).name()
            << endmsg;

    if( theDetectors.size()  != 1 ){

      std::stringstream es ;
      es << " getExtension: selection is not unique (or empty)  includeFlag: " << dd4hep::DetType( includeFlag ) << " excludeFlag: " << dd4hep::DetType( excludeFlag )
        << " --- found detectors : " ;
      for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ){
        es << theDetectors.at(i).name() << ", " ;
      }
      throw std::runtime_error( es.str() ) ;
    }

    theExtension = theDetectors.at(0).extension<dd4hep::rec::LayeredCalorimeterData>();

    return theExtension;
  }

  StatusCode initialize() override {

    // retrieve B field
    m_Bz = getFieldFromCompact();
    debug() << "B field (T) is : " << m_Bz << endmsg;

    // retrieve ecal dimensions:
    // - barrel: inner R, zmax
    // - endcap: inner R, zmin, zmax
    if (m_extrapolateToECal) {
      try {
        const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                                                                                        ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );
        m_eCalBarrelInnerR = eCalBarrelExtension->extent[0] / dd4hep::mm;
        m_eCalBarrelMaxZ = eCalBarrelExtension->extent[3] / dd4hep::mm;
        debug() << "ECAL barrel extent: Rmin [mm] = " << m_eCalBarrelInnerR << endmsg;
        debug() << "ECAL barrel extent: Zmax [mm] = " << m_eCalBarrelMaxZ << endmsg;
      }
      catch(...) {
        warning() << "ECAL barrel extension not found" << endmsg;
        m_eCalBarrelInnerR = 0.; // set to 0, will use it later to avoid projecting to the barrel
      };

      try {
        const dd4hep::rec::LayeredCalorimeterData * eCalEndCapExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                                                                                        ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );
        m_eCalEndCapInnerR = eCalEndCapExtension->extent[0] / dd4hep::mm;
        m_eCalEndCapOuterR = eCalEndCapExtension->extent[1] / dd4hep::mm;
        m_eCalEndCapInnerZ = eCalEndCapExtension->extent[2] / dd4hep::mm;
        m_eCalEndCapOuterZ = eCalEndCapExtension->extent[3] / dd4hep::mm;
        debug() << "ECAL endcap extent: Rmin [mm] = " << m_eCalEndCapInnerR << endmsg;
        debug() << "ECAL endcap extent: Rmax [mm] = " << m_eCalEndCapOuterR << endmsg;
        debug() << "ECAL endcap extent: Zmin [mm] = " << m_eCalEndCapInnerZ << endmsg;
        debug() << "ECAL endcap extent: Zmax [mm] = " << m_eCalEndCapOuterZ << endmsg;
      }
      catch(...) {
        warning() << "ECAL endcap extension not found" << endmsg;
        m_eCalEndCapInnerR = 0.; // set to 0, will use it later to avoid projecting to the endcap
      };
    }

    // setup system decoder
    m_systemEncoder = new dd4hep::DDSegmentation::BitFieldCoder(m_systemEncoding);
    m_indexSystem = m_systemEncoder->index("system");

    return StatusCode::SUCCESS;
  }

  std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection> operator()(const edm4hep::MCParticleCollection& genParticleColl, const SimTrackerHitColl& simTrackerHitCollVec) const override {

    // Create the output track collection and track gen<->reco links collection
    auto outputTrackCollection = edm4hep::TrackCollection();
    auto MCRecoTrackParticleAssociationCollection = edm4hep::TrackMCParticleLinkCollection();

    // loop over the gen particles, find charged ones, and create the corresponding reco particles
    int iparticle = 0;
    for (const auto& genParticle : genParticleColl) {
      debug() << endmsg;
      debug() << "Gen. particle: " << genParticle << endmsg;
      debug() <<"  particle "<<iparticle++<<"  PDG: "<< genParticle.getPDG()  << " energy: "<<genParticle.getEnergy()
              << " charge: "<< genParticle.getCharge() << endmsg;
      debug() << "Particle decayed in tracker: " << genParticle.isDecayedInTracker() << endmsg;

      // consider only charged particles
      if(genParticle.getCharge() == 0) continue;

      // GM: should skip also particles that cannot be reconstructed in tracker (vertex at too large radius, low energy)
      if (genParticle.getEnergy() < 0.010) continue; // cut at 10 MeV

      // Building an helix out of MCParticle properties and B field
      auto helixFromGenParticle = HelixClass_double();
      auto vertex = genParticle.getVertex();
      auto endpoint = genParticle.getEndpoint();
      debug() << "Vertex radius: " << sqrt(vertex.x*vertex.x+vertex.y*vertex.y) << endmsg;
      debug() << "Endpoint radius: " << sqrt(endpoint.x*endpoint.x+endpoint.y*endpoint.y) << endmsg;
      double genParticleVertex[] = {vertex.x, vertex.y, vertex.z};
      double genParticleMomentum[] = {genParticle.getMomentum().x, genParticle.getMomentum().y, genParticle.getMomentum().z};
      helixFromGenParticle.Initialize_VP(genParticleVertex, genParticleMomentum, genParticle.getCharge(), m_Bz);

      // Setting the track and trackStates at IP properties
      auto trackFromGen = edm4hep::MutableTrack();
      auto trackState_IP = edm4hep::TrackState {};
      trackState_IP.location = edm4hep::TrackState::AtIP;
      trackState_IP.D0 = helixFromGenParticle.getD0();
      trackState_IP.phi = helixFromGenParticle.getPhi0();
      trackState_IP.omega = helixFromGenParticle.getOmega();
      trackState_IP.Z0 = helixFromGenParticle.getZ0();
      trackState_IP.tanLambda = helixFromGenParticle.getTanLambda();
      trackState_IP.referencePoint = edm4hep::Vector3f((float)genParticleVertex[0],(float)genParticleVertex[1],(float)genParticleVertex[2]);
      // check this on CLD root files! Does it make any difference?
      // trackState_IP.referencePoint = edm4hep::Vector3f(0.0, 0.0, 0.0);
      trackFromGen.addToTrackStates(trackState_IP);

      // find SimTrackerHits associated to genParticle, store hit position, momentum and time
      // and calculate number of hits in each subdetector
      std::vector<std::array<double,7> > trackHits;
      std::vector<int32_t> v(m_trackerIDs.size());
      for ( size_t ih=0; ih<simTrackerHitCollVec.size(); ih++ ) {
        const edm4hep::SimTrackerHitCollection* coll = simTrackerHitCollVec[ih];
        for (const auto& hit : *coll) {
          const edm4hep::MCParticle particle = hit.getParticle();
          std::array<double,7> ahit{hit.x(), hit.y(), hit.z(), hit.getMomentum()[0], hit.getMomentum()[1], hit.getMomentum()[2], hit.getTime()};
          if(particle.getObjectID() == genParticle.getObjectID()) {
            trackHits.push_back(ahit);
            uint cellID = hit.getCellID();
            uint systemID = m_systemEncoder->get(cellID, m_indexSystem);
            for (size_t idxTracker=0; idxTracker < m_trackerIDs.size(); idxTracker++) {
              if (systemID == m_trackerIDs[idxTracker]) {
                v[idxTracker]++;
                break;
              }
            }
          }
        }
      }


      if(!trackHits.empty())
      {
        // particles with at least one SimTrackerHit
        debug() << "Number of SimTrackerHits: " << trackHits.size() << endmsg;

        // sort the hits according to their time
        std::sort(trackHits.begin(), trackHits.end(), [](const std::array<double,7> a, const std::array<double,7> b) {
          return a[6]<b[6];
        });

        // TrackState at First Hit
        auto trackState_AtFirstHit = edm4hep::TrackState {};
        auto firstHit = trackHits.front();
        double posAtFirstHit[] = {firstHit[0], firstHit[1], firstHit[2]};
        double momAtFirstHit[] = {firstHit[3], firstHit[4], firstHit[5]};
        debug() << "First hit: x, y, z, r = " << firstHit[0] << " "
                                              << firstHit[1] << " "
                                              << firstHit[2] << " "
                                              << sqrt(firstHit[0]*firstHit[0] + firstHit[1]*firstHit[1]) << endmsg;
        // get extrapolated momentum from the helix with ref point at IP
        helixFromGenParticle.getExtrapolatedMomentum(posAtFirstHit,momAtFirstHit);
        // produce new helix at first hit position
        auto helixAtFirstHit = HelixClass_double();
        helixAtFirstHit.Initialize_VP(posAtFirstHit, momAtFirstHit, genParticle.getCharge(), m_Bz);
        // fill the TrackState parameters
        trackState_AtFirstHit.location = edm4hep::TrackState::AtFirstHit;
        trackState_AtFirstHit.D0 = helixAtFirstHit.getD0();
        trackState_AtFirstHit.phi = helixAtFirstHit.getPhi0();
        trackState_AtFirstHit.omega = helixAtFirstHit.getOmega();
        trackState_AtFirstHit.Z0 = helixAtFirstHit.getZ0();
        trackState_AtFirstHit.tanLambda = helixAtFirstHit.getTanLambda();
        trackState_AtFirstHit.referencePoint = edm4hep::Vector3f((float)posAtFirstHit[0],(float)posAtFirstHit[1],(float)posAtFirstHit[2]);
        trackFromGen.addToTrackStates(trackState_AtFirstHit);

        // TrackState at Last Hit
        auto trackState_AtLastHit = edm4hep::TrackState{};
        auto lastHit = trackHits.back();
        double posAtLastHit[] = {lastHit[0], lastHit[1], lastHit[2]};
        double momAtLastHit[] = {lastHit[3], lastHit[4], lastHit[5]};
        debug() << "Last hit: x, y, z, r = " << lastHit[0] << " "
                                             << lastHit[1] << " "
                                             << lastHit[2] << " "
                                             << sqrt(lastHit[0]*lastHit[0] + lastHit[1]*lastHit[1]) << endmsg;
        // get extrapolated momentum from the helix with ref point at first hit
        helixAtFirstHit.getExtrapolatedMomentum(posAtLastHit, momAtLastHit);
        // produce new helix at last hit position
        auto helixAtLastHit = HelixClass_double();
        helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, genParticle.getCharge(), m_Bz);
        // fill the TrackState parameters
        trackState_AtLastHit.location = edm4hep::TrackState::AtLastHit;
        trackState_AtLastHit.D0 = helixAtLastHit.getD0();
        trackState_AtLastHit.phi = helixAtLastHit.getPhi0();
        trackState_AtLastHit.omega = helixAtLastHit.getOmega();
        trackState_AtLastHit.Z0 = helixAtLastHit.getZ0();
        trackState_AtLastHit.tanLambda = helixAtLastHit.getTanLambda();
        trackState_AtLastHit.referencePoint = edm4hep::Vector3f((float)posAtLastHit[0],
                                                                (float)posAtLastHit[1],
                                                                (float)posAtLastHit[2]);
        // attach the TrackState to the track
        trackFromGen.addToTrackStates(trackState_AtLastHit);

        // TrackState at Calorimeter
        if (m_eCalBarrelInnerR>0. || m_eCalEndCapInnerR>0.) {
          auto trackState_AtCalorimeter = edm4hep::TrackState{};
          pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);

          // create helix to project
          const pandora::Helix helix(trackState_IP.phi,
                                     trackState_IP.D0,
                                     trackState_IP.Z0,
                                     trackState_IP.omega,
                                     trackState_IP.tanLambda,
                                     m_Bz);
          const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
          const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);

          // First project to endcap
          float minGenericTime(std::numeric_limits<float>::max());
          if (m_eCalEndCapInnerR>0) {
            (void)helix.GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint,
                                    bestECalProjection, minGenericTime);
            // GM: if the radius of the point on the plane corresponding to the endcap inner face is
            // lower than the inner radius of the calorimeter, we might want to ignore it
            // for the moment let's keep it as it might be useful for debugging the reconstruction
          }

          // Then project to barrel surface(s), and keep projection with lower arrival time
          if (m_eCalBarrelInnerR>0) {
            pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
            float genericTime(std::numeric_limits<float>::max());
            const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));
            if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime)) {
              minGenericTime = genericTime;
              bestECalProjection = barrelProjection;
            }
            // GM: again, if the Z of the point on the cylinder of the barrel is beyond the
            // max/min z of the detector, we might want to ignore it - but let's keep it
            // for the moment as it might be useful for debugging the reconstruction
          }

          // get extrapolated position
          double posAtCalorimeter[] = {bestECalProjection.GetX(), bestECalProjection.GetY(), bestECalProjection.GetZ()};
          debug() << "Projection at calo: x, y, z, r = " << posAtCalorimeter[0] << " "
                                                         << posAtCalorimeter[1] << " "
                                                         << posAtCalorimeter[2] << " "
                                                         << sqrt(posAtCalorimeter[0]*posAtCalorimeter[0] + posAtCalorimeter[1]*posAtCalorimeter[1]) << endmsg;

          // get extrapolated momentum from the helix with ref point at last hit
          double momAtCalorimeter[] = {0.,0.,0.};
          helixAtLastHit.getExtrapolatedMomentum(posAtCalorimeter, momAtCalorimeter);

          // produce new helix at calorimeter position
          auto helixAtCalorimeter = HelixClass_double();
          helixAtCalorimeter.Initialize_VP(posAtCalorimeter, momAtCalorimeter, genParticle.getCharge(), m_Bz);

          // fill the TrackState parameters
          trackState_AtCalorimeter.location = edm4hep::TrackState::AtCalorimeter;
          trackState_AtCalorimeter.D0 = helixAtCalorimeter.getD0();
          trackState_AtCalorimeter.phi = helixAtCalorimeter.getPhi0();
          trackState_AtCalorimeter.omega = helixAtCalorimeter.getOmega();
          trackState_AtCalorimeter.Z0 = helixAtCalorimeter.getZ0();
          trackState_AtCalorimeter.tanLambda = helixAtCalorimeter.getTanLambda();
          trackState_AtCalorimeter.referencePoint = edm4hep::Vector3f((float)posAtCalorimeter[0],
                                                                      (float)posAtCalorimeter[1],
                                                                      (float)posAtCalorimeter[2]);
          // attach the TrackState to the track
          trackFromGen.addToTrackStates(trackState_AtCalorimeter);
        }

        // fill information about number of hits in the various subdetectors
        for (auto nhits : v) {
          trackFromGen.addToSubdetectorHitNumbers(nhits);
        }

        outputTrackCollection.push_back(trackFromGen);

        // Building the association between tracks and genParticles
        auto MCRecoTrackParticleAssociation = edm4hep::MutableTrackMCParticleLink();
        MCRecoTrackParticleAssociation.setFrom(trackFromGen);
        MCRecoTrackParticleAssociation.setTo(genParticle);
        MCRecoTrackParticleAssociationCollection.push_back(MCRecoTrackParticleAssociation);
      }
    }
    // push the output collections to event store
    return std::make_tuple(std::move(outputTrackCollection), std::move(MCRecoTrackParticleAssociationCollection));
  }

private:
  /// Solenoid magnetic field, to be retrieved from detector
  float m_Bz;

  /// ECAL barrel and endcap extent, to be retrieved from detector
  float m_eCalBarrelInnerR;
  float m_eCalBarrelMaxZ;
  float m_eCalEndCapInnerR;
  float m_eCalEndCapOuterR;
  float m_eCalEndCapInnerZ;
  float m_eCalEndCapOuterZ;

  /// configurable property to decide whether to calculate track state at ECAL or not
  Gaudi::Property<bool> m_extrapolateToECal{
    this, "ExtrapolateToECal", false, "Calculate track state at ECal inner face or not"
  };
  /// configurable properties (depend on detector) for calculating number of hits vs subdetector
  Gaudi::Property<std::string> m_systemEncoding{
    this, "SystemEncoding", "system:5", "System encoding string"
  };
  Gaudi::Property<std::vector<uint>> m_trackerIDs{
    this, "TrackerIDs", {}, "System IDs of tracking subdetectors"
  };

  /// General decoder to encode the tracker sub-system to determine
  dd4hep::DDSegmentation::BitFieldCoder* m_systemEncoder;
  int m_indexSystem;
};

DECLARE_COMPONENT(TracksFromGenParticles)
