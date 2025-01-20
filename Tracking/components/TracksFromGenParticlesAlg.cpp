// Gaudi
#include "Gaudi/Property.h"

// edm4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

// marlin
#include <marlinutil/HelixClass_double.h>

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// DD4HEP
#include "DD4hep/Detector.h"

// C++
#include <string>

/** @class TracksFromGenParticlesAlg
 *
 *  GaudiAlg version of TracksFromGenParticles, that builds an edm4hep::TrackCollection out of an edm4hep::MCParticleCollection.
 *  It just builds an helix out of the genParticle position, momentum, charge and user defined z component of the (constant) magnetic field.
 *  From this helix, different edm4hep::TrackStates (AtIP, AtFirstHit, AtLastHit and AtCalorimeter) are defined.
 *  The first and last hits are defined as those with smallest and largest R in the input SimTrackerHit collections
 *  This is meant to enable technical development needing edm4hep::Track and performance studies where having generator based trackis is a reasonable approximation.
 *  Possible inprovement:
 *    - Retrieve magnetic field from geometry: const DD4hep::Field::MagneticField& magneticField = detector.field(); DD4hep::DDRec::Vector3D field = magneticField.magneticField(point);
 *    - Handle properly tracks that hit the ECAL endcaps
 *  @author Brieuc Francois
 *  @author Archil Durglishvili
 *  @author Giovanni Marchiori
 */

class TracksFromGenParticlesAlg : public Gaudi::Algorithm {

  public:
    TracksFromGenParticlesAlg(const std::string& name, ISvcLocator* svcLoc);
    StatusCode initialize();
    StatusCode execute(const EventContext&) const;
    StatusCode finalize();

  private:
    /// Handle for input MC particles
    mutable DataHandle<edm4hep::MCParticleCollection> m_inputMCParticles{"MCParticles", Gaudi::DataHandle::Reader, this};
    /// List of input sim tracker hit collections
    Gaudi::Property<std::vector<std::string>> m_inputSimTrackerHitCollections{this, "InputSimTrackerHits", {}, "Names of SimTrackerHit collections to read"};
    /// the vector of input DataHandles for the tracker hit collections
    std::vector<DataObjectHandleBase*> m_inputSimTrackerHitCollectionHandles;
    /// Solenoid magnetic field, to be retrieved from detector
    float m_Bz;
    /// Inner radius of the ECAL
    Gaudi::Property<float> m_RadiusAtCalo{this, "RadiusAtCalo", 2172.8, "Inner radius (in mm) of the calorimeter where the TrackState::AtCalorimeter should be defined."};
    /// Handle for the output track collection
    mutable DataHandle<edm4hep::TrackCollection> m_tracks{"TracksFromGenParticlesAlg", Gaudi::DataHandle::Writer, this};
    /// Handle for the output links between reco and gen particles
    mutable DataHandle<edm4hep::TrackMCParticleLinkCollection> m_links{"TracksFromGenParticlesAlgAssociation", Gaudi::DataHandle::Writer, this};
};

TracksFromGenParticlesAlg::TracksFromGenParticlesAlg(const std::string& name, ISvcLocator* svcLoc) :
Gaudi::Algorithm(name, svcLoc) {
  declareProperty("InputGenParticles", m_inputMCParticles, "input MCParticles");
  declareProperty("OutputTracks", m_tracks, "Output tracks");
  declareProperty("OutputMCRecoTrackParticleAssociation", m_links, "MCRecoTrackParticleAssociation");

  m_Bz = 0.;
}

StatusCode TracksFromGenParticlesAlg::initialize() {
  StatusCode sc = Gaudi::Algorithm::initialize();
  if (sc.isFailure()) return sc;

  // FIXME: handle exceptions if collections not found
  for ( const auto& col : m_inputSimTrackerHitCollections ) {
    debug() << "Creating handle for input SimTrackerHit collection : " << col << endmsg;
    m_inputSimTrackerHitCollectionHandles.push_back(new DataHandle<edm4hep::SimTrackerHitCollection>(col, Gaudi::DataHandle::Reader, this));
  }

  // retrieve B field
  dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
  const double position[3]  = {0, 0, 0};  // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3] = {0, 0, 0};               // initialise object to hold magnetic field
  mainDetector.field().magneticField(position, magneticFieldVector);  // get the magnetic field vector from DD4hep
  m_Bz = magneticFieldVector[2] / dd4hep::tesla;  // z component at (0,0,0)
  debug() << "B field (T) is : " << m_Bz << endmsg;
  return StatusCode::SUCCESS;
}

StatusCode TracksFromGenParticlesAlg::execute(const EventContext&) const {
    // Get the input MC particle collection
    const edm4hep::MCParticleCollection* genParticleColl = m_inputMCParticles.get();
    // Create the output track collection and track gen<->reco links collection
    auto outputTrackCollection = new edm4hep::TrackCollection();
    auto MCRecoTrackParticleAssociationCollection = new edm4hep::TrackMCParticleLinkCollection();

    // loop over the gen particles, find charged ones, and create the corresponding reco particles
    int iparticle = 0;
    for (const auto& genParticle : *genParticleColl) {
      debug() << endmsg;
      debug() << "Gen. particle: " << genParticle << endmsg;
      debug() <<"  particle "<<iparticle++<<"  PDG: "<< genParticle.getPDG()  << " energy: "<<genParticle.getEnergy()
              << " charge: "<< genParticle.getCharge() << endmsg;
      debug() << "Particle decayed in tracker: " << genParticle.isDecayedInTracker() << endmsg;

      // consider only charged particles
      if (genParticle.getCharge() == 0) continue;

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
      trackFromGen.addToTrackStates(trackState_IP);

      // find SimTrackerHits associated to genParticle
      std::vector<std::array<double,6> > trackHits;
      for ( size_t ih=0; ih<m_inputSimTrackerHitCollectionHandles.size(); ih++ ) {
        auto handle = dynamic_cast<DataHandle<edm4hep::SimTrackerHitCollection>*> (m_inputSimTrackerHitCollectionHandles[ih]);
        const edm4hep::SimTrackerHitCollection* coll = handle->get();
        for (const auto& hit : *coll) {
          const edm4hep::MCParticle particle = hit.getParticle();
          std::array<double,6> ahit{hit.x(), hit.y(), hit.z(), hit.getMomentum()[0], hit.getMomentum()[1], hit.getMomentum()[2]};
          if(particle.getVertex() == genParticle.getVertex() && particle.getMomentum() == genParticle.getMomentum()) trackHits.push_back(ahit);
        }
      }
    
      if(!trackHits.empty())
      {
        // particles with at least one SimTrackerHit
        debug() << "Number of SimTrackerHits: " << trackHits.size() << endmsg;

        // sort the hits according to radial distance from beam axis
        // FIXME: does this work well also for tracks going to endcaps of vtx/wrapper?
        std::sort(trackHits.begin(), trackHits.end(), [](const std::array<double,6> a, const std::array<double,6> b) {
          double rhoA = std::sqrt(a[0]*a[0] + a[1]*a[1]);
          double rhoB = std::sqrt(b[0]*b[0] + b[1]*b[1]);
          return rhoA < rhoB;
        });

        // TrackState at First Hit
        auto trackState_AtFirstHit = edm4hep::TrackState {};
        auto firstHit = trackHits.front();
        double posAtFirstHit[] = {firstHit[0], firstHit[1], firstHit[2]};
        double momAtFirstHit[] = {firstHit[3], firstHit[4], firstHit[5]};
        debug() << "Radius of first hit: " << std::sqrt(firstHit[0]*firstHit[0] + firstHit[1]*firstHit[1]) << endmsg;
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
        debug() << "Radius of last hit: " << std::sqrt(lastHit[0]*lastHit[0] + lastHit[1]*lastHit[1]) << endmsg;
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
        auto trackState_AtCalorimeter = edm4hep::TrackState{};
        double pointAtCalorimeter[] = {0.,0.,0.,0.,0.,0.};
        auto time = helixAtLastHit.getPointOnCircle(m_RadiusAtCalo, posAtLastHit, pointAtCalorimeter);
        double posAtCalorimeter[] = {pointAtCalorimeter[0], pointAtCalorimeter[1], pointAtCalorimeter[2]};
        debug() << "Radius at calorimeter: " << std::sqrt(pointAtCalorimeter[0]*pointAtCalorimeter[0] + pointAtCalorimeter[1]*pointAtCalorimeter[1]) << endmsg;
        double momAtCalorimeter[] = {0.,0.,0.};
        // get extrapolated momentum from the helix with ref point at last hit
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

        outputTrackCollection->push_back(trackFromGen);

        // Building the association between tracks and genParticles
        auto MCRecoTrackParticleAssociation = edm4hep::MutableTrackMCParticleLink();
        MCRecoTrackParticleAssociation.setFrom(trackFromGen);
        MCRecoTrackParticleAssociation.setTo(genParticle);
        MCRecoTrackParticleAssociationCollection->push_back(MCRecoTrackParticleAssociation);
      }
    }

    // push the outputTrackCollection to event store
    m_tracks.put(outputTrackCollection);
    m_links.put(MCRecoTrackParticleAssociationCollection);

    debug() << "Output tracks collection size: " << outputTrackCollection->size() << endmsg;

    return StatusCode::SUCCESS;
}

StatusCode TracksFromGenParticlesAlg::finalize() { return Gaudi::Algorithm::finalize(); }


DECLARE_COMPONENT(TracksFromGenParticlesAlg)
