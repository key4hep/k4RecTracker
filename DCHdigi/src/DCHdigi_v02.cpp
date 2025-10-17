#include "DCHdigi_v02.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

// edm4hep
#include "edm4hep/utils/vector_utils.h"

// DD4hep
#include "DD4hep/Detector.h"

// STL
#include <random>
#include <ranges>



namespace {
    // A struct for saving the information required for calculating the number of clusters for each 
    // individual particle in one cell, acconting for different particle types producing different number of clusters
    struct ParticleClusterInfo {
        double beta_gamma = 0.0;
        double path_length_mm = 0.0;
    };

    struct HitInfo {
        const edm4hep::SimTrackerHit* simhit;
        double arrival_time_ns = 0.0;
        edm4hep::Vector3d position_mm;
        double distance_to_wire_mm = 0.0;
    };
}

DCHdigi_v02::DCHdigi_v02(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(
        name, 
        svcLoc,
        {
            KeyValues("InputSimHitCollection", {""}),
            KeyValues("HeaderName", {"EventHeader"})
        },
        {
            KeyValues("OutputDigihitCollection", {"DCHDigi2Collection"}),
            KeyValues("OutputLinkCollection", {"DCHDigi2SimLinkCollection"})
        }
    ) {}

StatusCode DCHdigi_v02::initialize() {

    m_uniqueIDSvc = serviceLocator()->service(m_uidSvcName);
    if (!m_uniqueIDSvc) {
        error() << "Unable to locate UniqueIDGenSvc with name: " << m_uidSvcName << endmsg;
        return StatusCode::FAILURE;
    }

    m_geoSvc = serviceLocator()->service(m_geoSvcName);
    if (!m_geoSvc) {
        error() << "Unable to locate GeoSvc with name: " << m_geoSvcName << endmsg;
        return StatusCode::FAILURE;
    }

    // Retrieve the subdetector
    std::string dch_name(m_dch_name.value());
    if (m_geoSvc->getDetector()->detectors().count(dch_name) == 0) {
        error() << "Detector <<" << dch_name << ">> does not exist." << endmsg;
        return StatusCode::FAILURE;
    }

    // Retrieve the detector element
    dd4hep::DetElement dch_detelem = m_geoSvc->getDetector()->detectors().at(dch_name);
    // Retrieve the DCH_info data extension for the drift chamber
    m_dch_info = dch_detelem.extension<dd4hep::rec::DCH_info>();
    if (not m_dch_info->IsValid()) {
        error() << "No valid data extension was found for detector <<" << dch_name << ">>." << endmsg;
        return StatusCode::FAILURE;
    }

    // Retrieve the readout associated with the detector element (subdetector)
    dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(dch_name);
    if (not dch_sd.isValid()) {
        error() << "No valid Sensitive Detector was found for detector <<" << dch_name << ">>." << endmsg;
        return StatusCode::FAILURE;
    }
    // set the cellID decoder
    m_decoder = dch_sd.readout().idSpec().decoder();

    m_event_counter = 0;
    return StatusCode::SUCCESS;
}

std::tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection> 
DCHdigi_v02::operator()(const edm4hep::SimTrackerHitCollection& input,
                        const edm4hep::EventHeaderCollection& header) const {

    extension::SenseWireHitCollection output;
    extension::SenseWireHitSimTrackerHitLinkCollection links;

    std::mt19937_64 random_engine;
    auto engine_seed = m_uniqueIDSvc->getUniqueID(header, this->name());
    random_engine.seed(engine_seed);

    // Create the random distributions for smearing the coordinates in dd4hep default units
    std::normal_distribution<double> gauss_z_ddu{0., m_z_resolution_mm.value() * dd4hep::mm};
    std::normal_distribution<double> gauss_xy_ddu{0., m_xy_resolution_mm.value() * dd4hep::mm};

    
    debug() << "Processing event " << ++m_event_counter << endmsg;
    debug() << "Processing SimTrackerHitCollection with " << input.size() << " hits." << endmsg;

    std::unordered_map<uint64_t, std::vector<edm4hep::SimTrackerHit>> cell_map;
    cell_map.reserve(input.size());
    // Loop over the inputs and save the cellIDs and corresponding hits to a map
    for (const auto& simhit: input){
        uint64_t cellID = simhit.getCellID();
        cell_map[cellID].push_back(simhit);
    }

    debug() << "Map contains " << cell_map.size() << " unique cellIDs." << endmsg;

    // Loop over the cells and do the digitisation of the simhits in each cell
    for (const auto& [cellID, simhits] : cell_map) {

        // Some geometry values needed for the calculations below
        int layer = m_dch_info->CalculateILayerFromCellIDFields(
            m_decoder->get(cellID, "layer"),
            m_decoder->get(cellID, "superlayer")
        );
        int nphi = m_decoder->get(cellID, "nphi");

/* THE FOLLOWING CALCULATION OF WIRE ANGLES HAS BEEN COPIED AS IS FROM DCHdigi_v01! */
        // The direction of the sense wires can be calculated as:
        //   RotationZ(WireAzimuthalAngle) * RotationX(stereoangle)
        // One point of the wire is for example the following:
        //   RotationZ(WireAzimuthalAngle) * Position(cell_rave_z0, 0 , 0)
        // variables aredefined below
        auto WireAzimuthalAngle = this->m_dch_info->Get_cell_phi_angle(layer, nphi);
        float WireStereoAngle = 0;
        {
            auto l = this->m_dch_info->database.at(layer);
            // radial middle point of the cell at Z=0
            auto cell_rave_z0 = 0.5 * (l.radius_fdw_z0 + l.radius_fuw_z0);
            // when building the twisted tube, the twist angle is defined as:
            //     cell_twistangle    = l.StereoSign() * DCH_i->twist_angle
            // which forces the stereoangle of the wire to have the oposite sign
            WireStereoAngle = (-1.) * l.StereoSign() * m_dch_info->stereoangle_z0(cell_rave_z0);
        }
/* END OF COPYING WIRE ANGLE CALCULATION FROM DCHdigi_v01 */

        // Hits will be grouped by time into batches separated by 400 ns (dead time of a cell)
        // Need an additional loop over the simhits to create these hit trains
        // Store relevant information for the hits in this cell in a vector
        std::vector<HitInfo> hit_info_vector;
        hit_info_vector.reserve(simhits.size());

        for (const auto& simhit : simhits) {
            //////////////////////////
            // POSITION INFORMATION //
            //////////////////////////

            // Get hit position to calculate distance to wire
            // Need to convert to TVector3 to use the DCH_info methods
            // Use dd4hep:mm as scale to convert into the dd4hep default units (_ddu)
            auto simhit_position_ddu = this->Convert_EDM4hepVector_to_TVector3(simhit.getPosition(), dd4hep::mm);

            auto hit_to_wire_vector_ddu = m_dch_info->Calculate_hitpos_to_wire_vector(layer, nphi, simhit_position_ddu);
            auto hit_projection_on_the_wire_ddu = simhit_position_ddu + hit_to_wire_vector_ddu;
            double distance_to_wire_mm = hit_to_wire_vector_ddu.Mag() / dd4hep::mm; // Explicitly cast to mm, no matter what the default unit is

            ////////////////////////////////
            // POSITION AND TIME SMEARING //
            ////////////////////////////////
            // xy smearing
            double smearing_xy_mm = gauss_xy_ddu(random_engine) / dd4hep::mm; 
            double digihit_distance_to_wire_mm = std::max(0.0, distance_to_wire_mm + smearing_xy_mm);
            double drift_time_ns = this->get_drift_time(digihit_distance_to_wire_mm);

            // z smearing
            double smearing_z_ddu = gauss_z_ddu(random_engine);
            TVector3 wire_direction_ez_ddu = (m_dch_info->Calculate_wire_vector_ez(layer, nphi)).Unit();
            hit_projection_on_the_wire_ddu += smearing_z_ddu * wire_direction_ez_ddu; // Need to multiply smearing_z_ddu with dd4hep::mm to cast into default units

            // Limit z to values inside the drift chamber.
            // NB: can lead to a peak at z=Lhalf and z=-Lhalf
            hit_projection_on_the_wire_ddu.SetZ(
                std::clamp(hit_projection_on_the_wire_ddu.Z(), -(this->m_dch_info->Lhalf), this->m_dch_info->Lhalf)
            );

            // Convert to edm4hep vector and cast to mm
            edm4hep::Vector3d digihit_position_mm = this->Convert_TVector3_to_EDM4hepVector(hit_projection_on_the_wire_ddu, 1.0 / dd4hep::mm);

            double distance_to_readout_mm = (this->m_dch_info->Lhalf/dd4hep::mm - std::abs(digihit_position_mm.z))/std::cos(WireStereoAngle);
            double travel_time_ns = this->get_signal_travel_time(distance_to_readout_mm);

            double arrival_time_ns = simhit.getTime() + drift_time_ns + travel_time_ns;

            // Create a HitInfo object and add it to the vector
            HitInfo hit_info;
            hit_info.simhit = &simhit;
            hit_info.arrival_time_ns = arrival_time_ns;
            hit_info.position_mm = digihit_position_mm;
            hit_info.distance_to_wire_mm = digihit_distance_to_wire_mm;
            hit_info_vector.push_back(hit_info);
        }

        // Sort the hit_info vector by time so that we can determine the trains
        std::ranges::sort(hit_info_vector, {}, &HitInfo::arrival_time_ns);

        
        std::vector<std::vector<HitInfo>> hit_train_vector;
        hit_train_vector.reserve(hit_info_vector.size());
        if (!hit_info_vector.empty()){
            // Start the first train with the very first hit
            hit_train_vector.push_back({ hit_info_vector.front() });

            for (size_t i = 1; i < hit_info_vector.size(); ++i) {
                const auto& current_hit = hit_info_vector[i];
                const auto& last_hit_in_train = hit_train_vector.back().back();

                // Check if the current hit is within the dead time of the last hit in the train
                if (current_hit.arrival_time_ns - last_hit_in_train.arrival_time_ns < m_deadtime_ns.value()) {
                    // If yes, add it to the current train
                    hit_train_vector.back().push_back(current_hit);
                } else {
                    // If no, start a new train
                    hit_train_vector.emplace_back(1, current_hit);
                }
            }
        }

        // Now loop over the time trains and create a digihit for each train
        // Each train is a collection of hits where there is no time gap between hits larger than the dead time
        for (const auto& hit_train : hit_train_vector) {
            if (hit_train.empty()) continue; // Shouldn't really occur, but just to be safe

            // Sum of all energy deposits in the train
            double edep_sum_GeV = 0.0;

            // Map to store the information for cluster calculation
            std::unordered_map<podio::ObjectID, ParticleClusterInfo> cluster_info_map;
            cluster_info_map.reserve(hit_train.size());

            // Loop over the hits in the train to do:
            // 1) Sum the energy deposits
            // 2) Collect information for cluster calculation
            for (const auto& hit_info: hit_train) {
                const auto& simhit = hit_info.simhit;

                // Integrate the deposited energy
                edep_sum_GeV += simhit->getEDep();

                /////////////////////////
                // CLUSTER INFORMATION //
                /////////////////////////

                // Each particle will create a different number of clusters (different beta*gamma values)
                // So store the MCparticles in a map for later cluster calculations
                auto mcparticle = simhit->getParticle();                
                auto object_id = mcparticle.getObjectID();

                // Check if the particle is already in the map
                auto [it, new_entry] = cluster_info_map.try_emplace(object_id, ParticleClusterInfo{});
                auto& cluster_info = it->second;

                // NB: at the moment, hits from secondary particles (isProducedBySecondary() flag) cannot be treated perfectly accurate,
                // because the beta*gamma of that secondary particle cannot be retrieved
                // This will need to be fixed in DCHdigi_v03, with the full waveform digitisation

                if (new_entry) {
                    // As simplfication: just take the betagamma value from the first hit we encounter (unlikely to change much over the course of one cell)
                    double mass_GeV = mcparticle.getMass();
                    double momentum_GeV = edm4hep::utils::magnitude(mcparticle.getMomentum());
                    cluster_info.beta_gamma = momentum_GeV / mass_GeV;
                    cluster_info.path_length_mm = simhit->getPathLength();
                } else {
                    // If the particle is already present, update the existing entry
                    cluster_info.path_length_mm += simhit->getPathLength();
                }

            } // end of loop over hit_info vector in this entry of the hit_train

            /////////////////////////
            // CLUSTER CALCULATION //
            /////////////////////////

            // Total number of clusters in cell built from each particle's contribution
            unsigned int total_nclusters = 0;

            // Calculate the number of clusters for each particle in the cell
            for (const auto& [object_id, particle_info] : cluster_info_map) {
                // Get number of clusters per length from delphes
                // Output from delphes function is in 1/m, so to convert to 1/mm we need to scale accordingly
                double nclusters_per_mm = m_delphesTrkUtil.Nclusters(particle_info.beta_gamma, m_GasSel.value()) / 1000.0;
                double nclusters_mean = nclusters_per_mm * particle_info.path_length_mm;

                std::poisson_distribution<int> poisson_dist(nclusters_mean);
                total_nclusters += poisson_dist(random_engine);
            }


            /////////////////////////////
            // SAVING THE SENSEWIREHIT //
            /////////////////////////////
            
            // Use the first hit in the train for the time, position, and distance to wire
            const auto& first_hit = hit_train.front();
            auto sense_wire_hit = output.create();
            sense_wire_hit.setCellID(cellID);
            sense_wire_hit.setType(0);
            sense_wire_hit.setQuality(0);
            sense_wire_hit.setTime(first_hit.arrival_time_ns);
            sense_wire_hit.setEDep(edep_sum_GeV);
            sense_wire_hit.setEDepError(0.0);
            sense_wire_hit.setPosition(first_hit.position_mm);
            sense_wire_hit.setPositionAlongWireError(m_z_resolution_mm);
            sense_wire_hit.setWireAzimuthalAngle(WireAzimuthalAngle);
            sense_wire_hit.setWireStereoAngle(WireStereoAngle);
            sense_wire_hit.setDistanceToWire(first_hit.distance_to_wire_mm);
            sense_wire_hit.setDistanceToWireError(m_xy_resolution_mm);

            // Clusters are added to the SenseWireHit as vector member containing the number of electrons in each cluster
            // The length of this vector is the total number of clusters in the cell
            // For now, we do not calculate the size of each cluster, so we just fill it with a dummy value
            for (unsigned int i = 0; i < total_nclusters; ++i) {
                sense_wire_hit.addToNElectrons(-1);
            }

            /////////////////////
            // SAVING THE LINK //
            /////////////////////

            auto link = links.create();
            link.setFrom(sense_wire_hit);
            link.setTo(*(first_hit.simhit));

        } // end of loop over hit_train_vector
        
    } // end of loop over the cells

    return std::make_tuple(std::move(output), std::move(links));

}


double DCHdigi_v02::get_drift_time(double distance_to_wire_mm) const {
    // Calculate the drift time based on the distance to the wire
    // This is a preliminary implementation that needs to be updated with a realistic model
    // For now, we use a simple linear model with a constant drift velocity

    double drift_velocity_um_per_ns = 25.0;

    // Convert distance to wire from mm to um
    double distance_to_wire_um = distance_to_wire_mm * 1000.0;

    // Calculate the drift time in ns
    double drift_time_ns = distance_to_wire_um / drift_velocity_um_per_ns;
    return drift_time_ns;
}


double DCHdigi_v02::get_signal_travel_time(double distance_to_readout_mm) const {
    // Calculate the time it takes for the signal to travel along the wire to the readout electronics
    // Assume 2/3 of the speed of light for the signal propagation speed

    double speed_of_light_mm_per_ns = 299.792458;
    double signal_speed_mm_per_ns = speed_of_light_mm_per_ns * 2.0 / 3.0;

    // Calculate the signal travel time in ns
    double signal_travel_time_ns = distance_to_readout_mm / signal_speed_mm_per_ns;
    return signal_travel_time_ns;
}
