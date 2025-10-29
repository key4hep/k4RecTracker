/** @class DCHdigi_v02
 * Gaudi Algorithm for DCH digitization
 *
 * @author Andreas Loeschcke Centeno
 * @date   2025-10-09
 *
 * Gaudi MultiTransformer that digitises SimTrackerHits from a Drift Chamber to edm4hep::SenseWireHits (currently a custom hit type defined in DCHdigi/dataFormatExtension)
 * 
 * In comparison to DCHdigi_v01, this version will produce only one DigiHit per cell, combining all SimHits in the same cell (unless there is siginificant time difference between the SimHits, larger than the m_deadtime_ns parameter).
 * To do this, the hits in the cells are sorted by time, after adding a (simplified, to be updated in the future) drift time and time to reach the readout.
 * They are further separated in hit 'trains' if the time difference between two hits is larger than the deadtime of a cell.
 * Each hit train creates one DigiHit, with the time and the position coming from the first hit (sorted by time) in the train.
 * The energy deposit will be the sum of all hits in the train.
 * 
 * The functional also does a dN/dx calculation for the cell, based on the parametrisation in delphes.
 * This is done by summing all the step lengths in the cell and getting the beta*gamma of the particle and passing it to the delphes parametrisation.
 * In case of multiple particles in the same cell, the number of clusters is calculated for each particle and summed up.
 * Since this uses delphes functions directly (via delphes::TrkUtil class member), this functional has a dependency on delphes.
 * The number of electrons within one cluster is not calculated at the moment, and filled with a dummy value (-1).
 *
 * Note: Variables for quantities with units attached to them, have either the units stated explicitly in the name as suffix (e.g. _mm, _ns) or by giving the unit system in which they are in (e.g. dd4hep default units: _ddu)
 *
 * Inputs:
 *     - SimTrackerHitCollection (SimTrackerHits in the DCH)
 *     - EventHeaderCollection (for consistently seeding the random engine)
 * 
 * Properties:
 *     - @param m_dch_name The name of the drift chamber geometry, needed to get the decoder for the cellID
 *     - @param m_z_resolution_mm Spatial resolution in the direction along the wire, in mm
 *     - @param m_xy_resolution_mm Spatial resolution in the direction perpendicular to the wire, in mm
 *     - @param m_deadtime_ns Deadtime of a cell in ns, hit trains in the same cell separated by more than this time will create separate DigiHits
 *     - @param m_GasType Gas type for the delphes dN/dx calculation: 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar
 *     - @param m_ReadoutWindowStartTime_ns Start time of the readout window in ns
 *     - @param m_ReadoutWindowDuration_ns Duration of the readout window in ns
 *     - @param m_uidSvcName The name of the UniqueIDGenSvc instance, used to create seed for each event/run, ensuring reproducibility.
 *     - @param m_geoSvcName The name of the GeoSvc instance, needed to intialise the DCH_info class for geometry calculations
 *
 * Outputs:
 *     - SenseWireHitCollection: Digitised hits
 *     - SenseWireHitSimTrackerHitLinkCollection: Links between SenseWireHits and the first SimTrackerHit in the train 
 * 
 * LIMITATIONS: (status 09/10/2025)
 *     - The drift time calculation is preliminary and needs to be updated with a realistic model
 *     - Potentially also the time to reach the readout could be updated
 *     - Hits with the isProducedBySecondary flag are not treated accurately for the dNdx calculation
 *     - The number of electrons in each cluster is not calculated
 * 
 */

#pragma once

// Gaudi
#include <GaudiKernel/SmartIF.h>
#include <GaudiKernel/ISvcLocator.h>
#include "Gaudi/Accumulators.h"

// k4FWCore
#include "k4FWCore/Transformer.h"

// k4Interface
#include <k4Interface/IUniqueIDGenSvc.h>
#include <k4Interface/IGeoSvc.h>

// edm4hep
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"

// edm4hep extension
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"

// DD4hep
#include "DDSegmentation/BitFieldCoder.h"

// DDRec
#include "DDRec/DCH_info.h"

// delphes
#include "TrackCovariance/TrkUtil.h"

// ROOT
#include "TRandom3.h"

class DCHdigi_v02 final
    : public k4FWCore::MultiTransformer<std::tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection>
                                       (const edm4hep::SimTrackerHitCollection&, 
                                        const edm4hep::EventHeaderCollection&)> {

public:
    DCHdigi_v02(const std::string& name, ISvcLocator* svcLoc);

    std::tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection> 
        operator()(const edm4hep::SimTrackerHitCollection& input,
                   const edm4hep::EventHeaderCollection& header) const override;

    StatusCode initialize() override final;

    // Fast zero-truncated Poisson sampler
    inline int sample_zero_truncated_poisson(double lambda, TRandom3& gen) const {
        double u = gen.Uniform(std::exp(-lambda), 1.0);
        double t = -std::log(u);
        int k = gen.Poisson(lambda - t);
        return 1 + k;
    }

private:
    SmartIF<IUniqueIDGenSvc> m_uniqueIDSvc{nullptr};
    Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "UniqueIDGenSvc", "The name of the UniqueIDGenSvc instance"};

    SmartIF<IGeoSvc> m_geoSvc{nullptr};
    Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

    mutable Gaudi::Accumulators::Counter<Gaudi::Accumulators::atomicity::full, unsigned int> m_event_counter;

    TrkUtil m_delphesTrkUtil;

    dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

    // Detector name
    Gaudi::Property<std::string> m_dch_name{this, "DCH_name", "DCH_v2", "Name of the Drift Chamber detector"};

    // Drift chamber info extension for geometry calculations
    dd4hep::rec::DCH_info* m_dch_info{nullptr};

    // z resolution in mm
    Gaudi::Property<double> m_z_resolution_mm{
        this, 
        "zResolution_mm", 
        1.0,
        "Spatial resolution in the direction along the wire, in mm."};
    // xy resolution in mm
    Gaudi::Property<double> m_xy_resolution_mm{
        this, 
        "xyResolution_mm", 
        0.1,
        "Spatial resolution in the direction perpendicular to the wire, in mm."};

    // Deadtime of a cell in ns
    Gaudi::Property<double> m_deadtime_ns{
        this,
        "Deadtime_ns",
        400.0,
        "Deadtime of a cell in ns."
    };

    // Gas drift velocity in um/ns
    Gaudi::Property<double> m_drift_velocity_um_per_ns{
        this,
        "DriftVelocity_um_per_ns",
        -1.0,
        "Gas drift velocity in um/ns. If negative, automatically chosen based on GasType. Currently assumed constant for the drift time calculation."
    };

    // Signal velocity in the wire in mm/ns
    Gaudi::Property<double> m_signal_velocity_mm_per_ns{
        this,
        "SignalVelocity_mm_per_ns",
        TMath::C()*1e-6*2.0/3.0, // 2/3 of the speed of light in mm/ns (1e-6 is to convert from m/s to mm/ns)
        "Signal velocity in the wire in mm/ns. Default value: 2/3 of the speed of light."
    };

    // Gas mixture in the chamber
    Gaudi::Property<int> m_GasType{
        this,
        "GasType",
        {0},
        "Gas type: 0: He(90%)-Isobutane(10%), 1: pure He, 2: Ar(50%)-Ethane(50%), 3: pure Ar."};

    // Readout window start time in ns
    Gaudi::Property<double> m_ReadoutWindowStartTime_ns{
        this,
        "ReadoutWindowStartTime_ns",
        1.0,
        "Together with ReadoutWindowDuration_ns, defines the readout window. Any DigiHits with arrival time before ReadoutWindowStartTime_ns are discarded."
    };

    // Readout window duration in ns
    Gaudi::Property<double> m_ReadoutWindowDuration_ns{
        this,
        "ReadoutWindowDuration_ns",
        450.0,
        "Together with ReadoutWindowStartTime_ns, defines the readout window. Any DigiHits with arrival time after ReadoutWindowStartTime_ns + ReadoutWindowDuration_ns are discarded."
    };


    /// Convert EDM4hep Vector3d to TVector3
    TVector3 toTVector3(const edm4hep::Vector3d& v) const {
        return {v[0], v[1], v[2]};
    };
    /// Convert TVector3 to EDM4hep Vector3d
    edm4hep::Vector3d toEDM4hepVector(const TVector3& v) const {
        return {v.x(), v.y(), v.z()};
    };

    // /// Function to calculate the drift time from the distance to the wire
    double get_drift_time_ns(double distance_to_wire_mm) const;

    // /// Function to calculate the time it takes for the signal to travel from the wire to the readout electronics
    double get_signal_travel_time_ns(double distance_to_readout_mm) const;

    /// Function to get the drift velocity based on the gas type
    double get_default_drift_velocity_um_per_ns() const;

};


DECLARE_COMPONENT(DCHdigi_v02);
