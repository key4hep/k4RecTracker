#pragma once

// Gaudi
#include <GaudiKernel/SmartIF.h>
#include <GaudiKernel/ISvcLocator.h>

// k4Interface
#include <k4Interface/IUniqueIDGenSvc.h>
#include <k4Interface/IGeoSvc.h>

// k4FWCore
#include "k4FWCore/Transformer.h"

// edm4hep
#include "edm4hep/TrackMCParticleLinkCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/RecDqdxCollection.h"

// delphes
#include "TrackCovariance/TrkUtil.h"


/** @class TrackdNdxDelphesBased
 *
 *  Gaudi transformer that builds an edm4hep::RecDqdxCollection out of an edm4hep::TrackMCParticleLinkCollection.
 *  It uses the delphes parametrisation for deriving the dN/dx based on a reconstructed track:
 *  For a given gas mixture, there exist data points specifying the number of clusters per meter for a given beta*gamma of the particle.
 *  A spline is used to interpolate between these data points.
 *  Additionally, in delphes there is a function to calculate the track length inside a cylindrical volume.
 *  The number of clusters per meter and the track length are used to derive the average number of cluster expected for the track.
 *  The 'true' value is drawn from a Poissonian and stored in the RecDqdxCollection.
 *  This Gaudi functional is linked against delphes, i.e. the delphes functions are called directly via an instance of the delphes::TrkUtil class.
 *  The benefit of this approach is that any future updates in the delphes code are automatically available here.
 *  For limitations, see below.
 *  First status of this algorithm was presented here: https://indico.cern.ch/event/1556632/#72-parametrized-cluster-counti
 *  
 *  PLEASE NOTE :
 *  This algortihm is mostly meant as work-around to be able to provide dN/dx values for track in order to enable studies which require such information,
 *  while work on a more accurate dN/dx digitiser on a cell-basis (instead of track-basis) is being conducted.
 * 
 *  Inputs: 
 *      - TrackMCParticleCollection (tracks for the track length calculation, MCparticle for the number of cluster calculation)
 *      - EventHeaderCollection     (for consistently seeding the random engine)
 *  Properties:
 *      - NAMES of the parameters under which detector dimensions are stored in the XML file, required for the track length calculation
 *        (thus, the dimension are loaded dynamically, defaults currently to IDEA drift chamber parameter names)
 *      - Gas selection: Integer to select the gas mixture for which the cluster data is currently available
 *      - Fill factor: A factor to allow the dN/dx calculation for tracks which are not in a continous active tracking volume (e.g. straw tube tracker).
 *        The track length is scaled by this factor, thus reducing the number of clusters in the end by the same factor.
 *  Outputs:
 *      - RecDqdxCollection: To each track, a quantity dqdx is associated, currently (13.06.2025) with only a value (not using error and type flag)
 *                           If the calculation somehow fails (see limitations below), the dqdx value is set to a dummy value (-999)
 * 
 *  LIMITATIONS (Status 13.06.2025) :
 *      - Entirely dependent on the implementation of the delphes functions:
 *          - Cluster calculation only available for 4 gas mixtures
 *          - Cluster calculation only for beta*gamma values for a certain range.
 *              - This algorithm sets dN/dx to the dummy value for particles outside the good beta*gamma range (warning produced)
 *          - Track length calculation is not available for particles with a large tranvserse impact parameters
 *              - If the track length calculation fails, dN/dx set to dummy (with warning)
 *          - Track length calculation assumes a perfectly cylindrical drift volume
 *      - Since this parametrisation is based on full tracks, the energy loss in the tracking volume cannot be accounted for
 *      - dN/dx quantity is only filled with a value, no error at the moment
 *  
 *  @author Andreas Loeschcke Centeno
 */
class TrackdNdxDelphesBased final
    : public k4FWCore::Transformer<edm4hep::RecDqdxCollection(const edm4hep::TrackMCParticleLinkCollection&, const edm4hep::EventHeaderCollection&)> {

public:
    TrackdNdxDelphesBased(const std::string& name, ISvcLocator* svcLoc);

    edm4hep::RecDqdxCollection operator()(const edm4hep::TrackMCParticleLinkCollection& input,
                                          const edm4hep::EventHeaderCollection& header) const override;

    StatusCode initialize() final;

private:
    SmartIF<IUniqueIDGenSvc> m_uniqueIDSvc{nullptr};
    SmartIF<IGeoSvc> m_geoSvc{nullptr};

    TrkUtil m_delphesTrkUtil;

    Gaudi::Property<std::string> m_Zmax_parameter_name{
        this,
        "ZmaxParameterName",
        {"DCH_gas_Lhalf"},
        "Name of XML file parameter describing the z extent of the active volume (max value in +z direction)"};
    Gaudi::Property<std::string> m_Zmin_parameter_name{
        this,
        "ZminParameterName",
        {"DCH_gas_Lhalf"},
        "Name of XML file parameter describing the z extent of the active volume (min value in -z direction). For forward-backward symmetric detectors, use same name as for ZmaxParameterName. This value is then automatically converted into the negative value."};
    Gaudi::Property<std::string> m_Rmin_parameter_name{
        this,
        "RminParameterName",
        {"DCH_gas_inner_cyl_R"},
        "Name of XML file parameter describing the inner radius of the active volume."};
    Gaudi::Property<std::string> m_Rmax_parameter_name{
        this,
        "RmaxParameterName",
        {"DCH_gas_outer_cyl_R"},
        "Name of XML file parameter describing the outer radius of the active volume."};
    Gaudi::Property<int> m_GasSel{
        this,
        "GasSel",
        {0},
        "Gas selection: 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar."};
    Gaudi::Property<double> m_fill_factor{
        this,
        "FillFactor",
        {1.0},
        "Factor (between 0 and 1) describing the fraction of the detector volume that is active (e.g., for Straw Tube Tracker the factor is significantly below 1 due to gaps between tubes). The factor is used to scale the calculated track length and thus the corresponding number of clusters."};
};
