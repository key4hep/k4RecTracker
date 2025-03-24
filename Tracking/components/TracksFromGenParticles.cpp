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

  edm4hep::TrackState getExtrapolationAtCalorimeter(const pandora::CartesianVector& ecalProjection, const HelixClass_double& helixAtLastHit) const {

    edm4hep::TrackState trackState_AtCalorimeter = edm4hep::TrackState{};

    double posAtCalorimeter[] = {ecalProjection.GetX(), ecalProjection.GetY(), ecalProjection.GetZ()};
    debug()
      << "Projection at calo: x, y, z, r = "
      << posAtCalorimeter[0] << " "
      << posAtCalorimeter[1] << " "
      << posAtCalorimeter[2] << " "
      << sqrt(posAtCalorimeter[0]*posAtCalorimeter[0] + posAtCalorimeter[1]*posAtCalorimeter[1]) << endmsg;

    // get extrapolated momentum from the helix with ref point at last hit
    double momAtCalorimeter[] = {0.,0.,0.};
    helixAtLastHit.getExtrapolatedMomentum(posAtCalorimeter, momAtCalorimeter);

    // produce new helix at calorimeter position
    auto helixAtCalorimeter = HelixClass_double();
    helixAtCalorimeter.Initialize_VP(posAtCalorimeter, momAtCalorimeter, helixAtLastHit.getCharge(), m_Bz);

    // fill the TrackState parameters
    trackState_AtCalorimeter.location = edm4hep::TrackState::AtCalorimeter;
    trackState_AtCalorimeter.D0 = helixAtCalorimeter.getD0();
    trackState_AtCalorimeter.phi = std::atan2(momAtCalorimeter[1], momAtCalorimeter[0]);
    trackState_AtCalorimeter.omega = helixAtCalorimeter.getOmega();
    trackState_AtCalorimeter.Z0 = helixAtCalorimeter.getZ0();
    trackState_AtCalorimeter.tanLambda = helixAtCalorimeter.getTanLambda();
    trackState_AtCalorimeter.referencePoint = edm4hep::Vector3f(posAtCalorimeter[0],
                                                                posAtCalorimeter[1],
                                                                posAtCalorimeter[2]);
    return trackState_AtCalorimeter;
  }


  std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection> operator()(const edm4hep::MCParticleCollection& genParticleColl, const SimTrackerHitColl& simTrackerHitCollVec) const override {

    // Create the output track collection and track gen<->reco links collection
    auto outputTrackCollection = edm4hep::TrackCollection();
    auto MCRecoTrackParticleAssociationCollection = edm4hep::TrackMCParticleLinkCollection();

    // loop over the gen particles, find charged ones, and create the corresponding reco particles
    int iparticle = 0;
    for (const auto& genParticle : genParticleColl) {
      edm4hep::Vector3d p = genParticle.getMomentum();
      double pmag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
      debug() << endmsg;
      debug() << "Gen. particle: " << genParticle << endmsg;
      debug() <<"  particle "<<iparticle++<<"  PDG: "<< genParticle.getPDG()  << " momentum: "<< pmag
              << " charge: "<< genParticle.getCharge() << endmsg;
      debug() << "Particle decayed in tracker: " << genParticle.isDecayedInTracker() << endmsg;

      // consider only charged particles
      if(genParticle.getCharge() == 0) continue;

      // skip low momentum particles that cannot be reconstructed in tracker
      if (pmag < m_minParticleMomentum) continue;

      // build an helix out of MCParticle properties and B field
      auto helixFromGenParticle = HelixClass_double();
      auto vertex = genParticle.getVertex();
      auto endpoint = genParticle.getEndpoint();
      debug() << "Vertex radius: " << sqrt(vertex.x*vertex.x+vertex.y*vertex.y) << endmsg;
      debug() << "Endpoint radius: " << sqrt(endpoint.x*endpoint.x+endpoint.y*endpoint.y) << endmsg;
      double genParticleVertex[] = {vertex.x, vertex.y, vertex.z};
      double genParticleMomentum[] = {genParticle.getMomentum().x, genParticle.getMomentum().y, genParticle.getMomentum().z};
      helixFromGenParticle.Initialize_VP(genParticleVertex, genParticleMomentum, genParticle.getCharge(), m_Bz);

      // set the track and trackStates at IP properties
      auto trackFromGen = edm4hep::MutableTrack();
      auto trackState_IP = edm4hep::TrackState {};
      trackState_IP.location = edm4hep::TrackState::AtIP;
      trackState_IP.D0 = helixFromGenParticle.getD0();
      // note: the phi in the EDM track state is the phi of the momentum vector at reference, see
      // https://github.com/iLCSoft/DDMarlinPandora/blob/master/src/DDTrackCreatorBase.cc#L360
      trackState_IP.phi = std::atan2(genParticle.getMomentum().y, genParticle.getMomentum().x);
      trackState_IP.omega = helixFromGenParticle.getOmega();
      trackState_IP.Z0 = helixFromGenParticle.getZ0();
      trackState_IP.tanLambda = helixFromGenParticle.getTanLambda();
      trackState_IP.referencePoint = edm4hep::Vector3f(genParticleVertex[0], genParticleVertex[1], genParticleVertex[2]);
      trackFromGen.addToTrackStates(trackState_IP);

      // find SimTrackerHits associated to genParticle (and not produced by secondaries)
      // store hit position, momentum and time
      // and calculate number of hits in each subdetector
      std::vector<std::array<double,7> > trackHits;
      std::vector<int> v(m_trackerIDs.size());
      for ( size_t ih=0; ih<simTrackerHitCollVec.size(); ih++ ) {
        const edm4hep::SimTrackerHitCollection* coll = simTrackerHitCollVec[ih];
        for (const auto& hit : *coll) {
          // skip hits that are produced by secondary particles:
          // they are pointing to the parent particles if the secondary one
          // is not kept in the MCParticle collection
          if (hit.isProducedBySecondary()) continue;

          // check that the ID of the particle that created the it is the same as the MCParticle being considered
          const edm4hep::MCParticle particle = hit.getParticle();
          if(particle.getObjectID() == genParticle.getObjectID()) {

            // store hit position, track momentum at hit and hit time in trackHits
            std::array<double,7> ahit{hit.x(), hit.y(), hit.z(), hit.getMomentum()[0], hit.getMomentum()[1], hit.getMomentum()[2], hit.getTime()};
            trackHits.push_back(ahit);

            // find systemID of hit and increase hit counter for corresponding subdetector
            uint cellID = hit.getCellID();
            int systemID = m_systemEncoder->get(cellID, m_indexSystem);
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
        double rAtFirstHit = sqrt(firstHit[0]*firstHit[0] + firstHit[1]*firstHit[1]);
        debug() << "First hit: x, y, z, r = "
                << firstHit[0] << " "
                << firstHit[1] << " "
                << firstHit[2] << " "
                << rAtFirstHit << endmsg;
        debug() << "First hit: px, py, pz = "
                << firstHit[3] << " "
                << firstHit[4] << " "
                << firstHit[5] << endmsg;

        // produce new helix at first hit position
        auto helixAtFirstHit = HelixClass_double();
        helixAtFirstHit.Initialize_VP(posAtFirstHit, momAtFirstHit, genParticle.getCharge(), m_Bz);
        // fill the TrackState parameters
        trackState_AtFirstHit.location = edm4hep::TrackState::AtFirstHit;
        trackState_AtFirstHit.D0 = helixAtFirstHit.getD0();
        trackState_AtFirstHit.phi = std::atan2(momAtFirstHit[1], momAtFirstHit[0]);
        trackState_AtFirstHit.omega = helixAtFirstHit.getOmega();
        trackState_AtFirstHit.Z0 = helixAtFirstHit.getZ0();
        trackState_AtFirstHit.tanLambda = helixAtFirstHit.getTanLambda();
        trackState_AtFirstHit.referencePoint = edm4hep::Vector3f(posAtFirstHit[0], posAtFirstHit[1], posAtFirstHit[2]);
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
        debug() << "Last hit: px, py, pz = "
                << lastHit[3] << " "
                << lastHit[4] << " "
                << lastHit[5] << endmsg;

        // produce new helix at last hit position
        auto helixAtLastHit = HelixClass_double();
        helixAtLastHit.Initialize_VP(posAtLastHit, momAtLastHit, genParticle.getCharge(), m_Bz);
        // fill the TrackState parameters
        trackState_AtLastHit.location = edm4hep::TrackState::AtLastHit;
        trackState_AtLastHit.D0 = helixAtLastHit.getD0();
        trackState_AtLastHit.phi = std::atan2(momAtLastHit[1], momAtLastHit[0]);
        trackState_AtLastHit.omega = helixAtLastHit.getOmega();
        trackState_AtLastHit.Z0 = helixAtLastHit.getZ0();
        trackState_AtLastHit.tanLambda = helixAtLastHit.getTanLambda();
        trackState_AtLastHit.referencePoint = edm4hep::Vector3f(posAtLastHit[0], posAtLastHit[1], posAtLastHit[2]);
        // attach the TrackState to the track
        trackFromGen.addToTrackStates(trackState_AtLastHit);

        // TrackState at Calorimeter
        if (m_eCalBarrelInnerR>0. || m_eCalEndCapInnerR>0.) {
          pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
          pandora::CartesianVector secondBestECalProjection(0.f, 0.f, 0.f);
          float minGenericTime(std::numeric_limits<float>::max());

          // create helix to project
          // rather than using parameters at production, better to use those from
          // last hit
          pandora::CartesianVector pos(posAtLastHit[0], posAtLastHit[1], posAtLastHit[2]);
          pandora::CartesianVector mom(momAtLastHit[0], momAtLastHit[1], momAtLastHit[2]);
          const pandora::Helix helix(pos, mom, genParticle.getCharge(), m_Bz);
          const pandora::CartesianVector& referencePoint(helix.GetReferencePoint());
          const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);

          // First project to endcap
          pandora::CartesianVector endCapProjection(0.f, 0.f, 0.f);
          bool hasEndCapProjection(false);
          if (m_eCalEndCapInnerR>0) {
            float genericTime(std::numeric_limits<float>::max());
            const pandora::StatusCode statusCode(helix.GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint,
                                                 endCapProjection, genericTime));
            float x = endCapProjection.GetX();
            float y = endCapProjection.GetY();
            float r = std::sqrt(x*x+y*y);
            if (
              (pandora::STATUS_CODE_SUCCESS == statusCode) &&
              (genericTime < minGenericTime) &&
              (r >= m_eCalEndCapInnerR) &&
              (r <= m_eCalEndCapOuterR)
            ) {
              minGenericTime = genericTime;
              bestECalProjection = endCapProjection;
              hasEndCapProjection = true;
            }
          }

          // Then project to barrel surface(s), and keep projection
          // if extrapolation is within the z acceptance of the detector
          pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
          bool hasBarrelProjection = false;
          if (m_eCalBarrelInnerR>0) {
            float genericTime(std::numeric_limits<float>::max());
            const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_eCalBarrelInnerR, referencePoint,
                                                 barrelProjection, genericTime));
            if (
              (pandora::STATUS_CODE_SUCCESS == statusCode) &&
              (std::fabs(barrelProjection.GetZ())<= m_eCalBarrelMaxZ)
            ) {
              hasBarrelProjection = true;
              if (genericTime < minGenericTime) {
                minGenericTime = genericTime;
                secondBestECalProjection = bestECalProjection;
                bestECalProjection = barrelProjection;
              }
              else {
                secondBestECalProjection = barrelProjection;
              }
            }
          }

          // store extrapolation to calo
          // by default, store extrapolation with lower arrival time
          // get extrapolated position
          edm4hep::TrackState trackState_AtCalorimeter = getExtrapolationAtCalorimeter(bestECalProjection, helixAtLastHit);

          // attach the TrackState to the track
          trackFromGen.addToTrackStates(trackState_AtCalorimeter);

          // attach second extrapolation if desired
          if (!m_keepOnlyBestExtrapolation and hasBarrelProjection and hasEndCapProjection) {
            edm4hep::TrackState trackState_AtCalorimeter_2 = getExtrapolationAtCalorimeter(secondBestECalProjection, helixAtLastHit);
            trackState_AtCalorimeter_2.location = edm4hep::TrackState::AtOther;
            trackFromGen.addToTrackStates(trackState_AtCalorimeter_2);
          }
        }

        // fill information about number of hits in the various subdetectors
        for (auto nhits : v) {
          trackFromGen.addToSubdetectorHitNumbers(nhits);
        }

        // add track to output collection
        outputTrackCollection.push_back(trackFromGen);

        // build the association between tracks and genParticles
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

  /// Configurable property to decide whether to calculate track state at ECAL or not
  Gaudi::Property<bool> m_extrapolateToECal{
    this, "ExtrapolateToECal", false, "Calculate track state at ECal inner face or not"
  };

  /// Configurable property to keep only first extrapolation to ECAL or also second one
  /// if both barrel and endcap are reached. The second projection will be set with
  /// location AtOther as there reconstruction (LCIO preprocessing for Pandora) forbids
  /// a second AtCalorimeter track state
  Gaudi::Property<bool> m_keepOnlyBestExtrapolation{
    this, "KeepOnlyBestExtrapolation", true, "Keep only extrapolation with shortest time or not"
  };

  /// Configurable property to keep only particles with energy above a given threshold
  Gaudi::Property<float> m_minParticleMomentum{
    this, "MinimumParticleMomentum", 0.010, "Keep only particles with momentum (in GeV) greater than MinimumParticleMomentum"
  };

  /// Configurable property listing the systemIDs of the variuos tracker subdetectors
  Gaudi::Property<std::vector<int>> m_trackerIDs{
    this, "TrackerIDs", {}, "System IDs of tracking subdetectors"
  };

  /// General decoder to retrieve from each hit what is the
  /// system it belongs to. The tool will count number of hits in the different
  /// tracking subsystems based on the hit systemID, and on the list of systemIDs
  /// passed through TrackerIDs
  dd4hep::DDSegmentation::BitFieldCoder* m_systemEncoder;

  /// Configurable property storing string and number of bits used
  /// to encode the systemID in the hits of the various tracking devices
  /// (we do not care about the other fields of the readout)
  Gaudi::Property<std::string> m_systemEncoding{
    this, "SystemEncoding", "system:5", "System encoding string"
  };

  /// Used to retrieve systemID by index rather than by string
  int m_indexSystem;
};

DECLARE_COMPONENT(TracksFromGenParticles)
