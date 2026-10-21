#ifndef ELECTRON_DECORATOR_H
#define ELECTRON_DECORATOR_H

#include <k4FWCore/Transformer.h>

#include <edm4hep/ReconstructedParticleCollection.h>

// marlin
#include <marlinutil/HelixClass_double.h>

// pandora
#include <Objects/Helix.h>

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Volumes.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4Interface/IGeoSvc.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/TrackCollection.h"
#include <podio/LinkCollection.h>

class ElectronDecorator:
    public k4FWCore::Transformer<
        edm4hep::ReconstructedParticleCollection(
            const edm4hep::ReconstructedParticleCollection&)> {

public:

    ElectronDecorator(const std::string& name,
                    ISvcLocator* svc);

    edm4hep::ReconstructedParticleCollection operator()(
        const edm4hep::ReconstructedParticleCollection&) const override;

private:
     ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Geometry service"};

     Gaudi::Property<double> m_dEtaMargin{this, "dEtaMargin", 0.05, "Eta matching margin"};
     Gaudi::Property<double> m_dPhiMargin{this, "dPhiMargin", 0.005, "Phi matching margin"};
};

#endif // ELECTRON_DECORATOR_H
