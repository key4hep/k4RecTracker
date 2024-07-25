#include "Gaudi/Property.h"

// edm4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCRecoTrackParticleAssociationCollection.h"

// marlin
#include <marlinutil/HelixClass_double.h>

// k4FWCore
#include "k4FWCore/Transformer.h"

#include <string>

/** @class TracksFromGenParticles
 *
 *  Gaudi transformer that builds an edm4hep::TrackCollection out of an edm4hep::MCParticleCollection.
 *  It just builds an helix out of the genParticle position, momentum, charge and user defined z component of the (constant) magnetic field.
 *  From this helix, different edm4hep::TrackStates (AtIP, AtFirstHit, AtLastHit and AtCalorimeter) are defined. #FIXME for now these trackstates are dummy (copies of the same helix parameters)
 *  This is meant to enable technical development needing edm4hep::Track and performance studies where having generator based trackis is a reasonable approximation.
 *  Possible inprovement:
 *    - Retrieve magnetic field from geometry: const DD4hep::Field::MagneticField& magneticField = detector.field(); DD4hep::DDRec::Vector3D field = magneticField.magneticField(point);
 *    - Properly define different trackStates
 *
 *  @author Brieuc Francois
 */

struct TracksFromGenParticles final
  : Gaudi::Functional::MultiTransformer<std::tuple<edm4hep::TrackCollection, edm4hep::MCRecoTrackParticleAssociationCollection>(const edm4hep::MCParticleCollection&)> {
  TracksFromGenParticles(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(
            name, svcLoc,
            {KeyValues("InputGenParticles", {"MCParticles"})},
            {KeyValues("OutputTracks", {"TracksFromGenParticles"}),
            KeyValues("OutputMCRecoTrackParticleAssociation", {"TracksFromGenParticlesAssociation"})}) {
  }

std::tuple<edm4hep::TrackCollection, edm4hep::MCRecoTrackParticleAssociationCollection> operator()(const edm4hep::MCParticleCollection& genParticleColl) const override {

    auto outputTrackCollection = edm4hep::TrackCollection();
    auto MCRecoTrackParticleAssociationCollection = edm4hep::MCRecoTrackParticleAssociationCollection();

    for (const auto& genParticle : genParticleColl) {
      debug() << "Particle decayed in tracker: " << genParticle.isDecayedInTracker() << endmsg;
      debug() << genParticle << endmsg;

      // Building an helix out of MCParticle properties and B field
      auto helixFromGenParticle = HelixClass_double();
      double genParticleVertex[] = {genParticle.getVertex().x, genParticle.getVertex().y, genParticle.getVertex().z};
      double genParticleMomentum[] = {genParticle.getMomentum().x, genParticle.getMomentum().y, genParticle.getMomentum().z};
      helixFromGenParticle.Initialize_VP(genParticleVertex, genParticleMomentum, genParticle.getCharge(), m_Bz);

      // Setting the track and trackStates properties
      // #FIXME for now, the different trackStates are dummy
      auto trackFromGen = edm4hep::MutableTrack();
      auto trackState_IP = edm4hep::TrackState {};
      trackState_IP.location = edm4hep::TrackState::AtIP;
      trackState_IP.D0 = helixFromGenParticle.getD0();
      trackState_IP.phi = helixFromGenParticle.getPhi0();
      trackState_IP.omega = helixFromGenParticle.getOmega();
      trackState_IP.Z0 = helixFromGenParticle.getZ0();
      trackState_IP.tanLambda = helixFromGenParticle.getTanLambda();
      trackFromGen.addToTrackStates(trackState_IP);
      auto trackState_AtFirstHit = edm4hep::TrackState(trackState_IP);
      trackState_AtFirstHit.location = edm4hep::TrackState::AtFirstHit;
      trackFromGen.addToTrackStates(trackState_AtFirstHit);
      auto trackState_AtLastHit = edm4hep::TrackState(trackState_IP);
      trackState_AtFirstHit.location = edm4hep::TrackState::AtLastHit;
      trackFromGen.addToTrackStates(trackState_AtLastHit);
      auto trackState_AtCalorimeter = edm4hep::TrackState(trackState_IP);
      trackState_AtFirstHit.location = edm4hep::TrackState::AtCalorimeter;
      trackFromGen.addToTrackStates(trackState_AtCalorimeter);

      debug() << trackFromGen << endmsg;
      outputTrackCollection.push_back(trackFromGen);

      // Building the association between tracks and genParticles
      auto MCRecoTrackParticleAssociation = edm4hep::MutableMCRecoTrackParticleAssociation();
      MCRecoTrackParticleAssociation.setRec(trackFromGen);
      MCRecoTrackParticleAssociation.setSim(genParticle);
      MCRecoTrackParticleAssociationCollection.push_back(MCRecoTrackParticleAssociation);
    }
    return std::make_tuple(std::move(outputTrackCollection), std::move(MCRecoTrackParticleAssociationCollection));
  }

  Gaudi::Property<float> m_Bz{this, "Bz", 2., "Z component of the (assumed constant) magnetic field in Tesla."};
};

DECLARE_COMPONENT(TracksFromGenParticles)
