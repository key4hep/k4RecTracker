#include "DCHdigi_v02.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

// edm4hep
#include "edm4hep/utils/vector_utils.h"

// DD4hep
#include "DD4hep/Detector.h"

// STL
#include <random>
// #include <bitset>


// A struct for saving the information required for calculating the number of clusters for each 
// individual particle in one cell, acconting for different particle types producing different number of clusters
namespace {
    struct ParticleClusterInfo {
        double beta_gamma = 0.0;
        double path_length = 0.0;
    };
}

DCHdigi_v02::DCHdigi_v02(const std::string& name, ISvcLocator* svcLoc)
    : Transformer(
        name, 
        svcLoc,
        {
            KeyValues("InputSimHitCollection", {""}),
            KeyValues("HeaderName", {"EventHeader"})
        },
        {
            KeyValues("OutputCollection", {"DCHDigi2Collection"})
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

extension::SenseWireHitCollection DCHdigi_v02::operator()(const edm4hep::SimTrackerHitCollection& input,
                                                        const edm4hep::EventHeaderCollection& header) const {

    extension::SenseWireHitCollection output;

    std::mt19937_64 random_engine;
    auto engine_seed = m_uniqueIDSvc->getUniqueID(header, this->name());
    random_engine.seed(engine_seed);

    // Create the random distributions for smearing the coordinates in dd4hep default units
    std::normal_distribution<double> gauss_z{0., m_z_resolution.value() * dd4hep::mm};
    std::normal_distribution<double> gauss_xy{0., m_xy_resolution.value() * dd4hep::mm};

    
    debug() << "Processing event " << ++m_event_counter << endmsg;
    debug() << "Processing SimTrackerHitCollection with " << input.size() << " hits." << endmsg;
    std::unordered_map<uint64_t, std::vector<edm4hep::SimTrackerHit>> cell_map;
    // Loop over the inputs and save the cellIDs and corresponding hits to a map
    for (const auto& simhit: input){
        uint64_t cellID = simhit.getCellID();
        cell_map[cellID].push_back(simhit);
    }

    debug() << "Map contains " << cell_map.size() << " unique cellIDs." << endmsg;

    for (const auto& [cellID, simhits] : cell_map) {
        auto sense_wire_hit = output.create();

        // Save the closest hit to the wire
        // Digitised hit time and position will be based on closest hit
        const edm4hep::SimTrackerHit* closest_simhit = nullptr;
        // Save also the vectors that were calculated anyway to find the closest hit
        // They will be reused for digitisation of the hit
        TVector3 closest_simhit_position;
        TVector3 closest_hit_to_wire_vector;
        double closest_distance_to_wire_mm = std::numeric_limits<double>::max();
        double edep_sum = 0.0;

        // Finding hit closest to the wire requires global layer number, nphi
        int layer = m_dch_info->CalculateILayerFromCellIDFields(
            m_decoder->get(cellID, "layer"),
            m_decoder->get(cellID, "superlayer")
        );
        int nphi = m_decoder->get(cellID, "nphi");

        // Map to store the information for cluster calculation
        std::unordered_map<podio::ObjectID, ParticleClusterInfo> cluster_info_map;

        // Loop over the simhits in the cell
        for (const auto& simhit : simhits) {
            
            //////////////////////////
            // POSITION INFORMATION //
            //////////////////////////

            // Get hit position to calculate distance to wire
            // Need to convert to TVector3 to use the DCH_info methods
            // Use dd4hep:mm as scale to convert into the dd4hep default unit system
            auto simhit_position = this->Convert_EDM4hepVector_to_TVector3(simhit.getPosition(), dd4hep::mm);

            auto hit_to_wire_vector = m_dch_info->Calculate_hitpos_to_wire_vector(layer, nphi, simhit_position);
            double distance_to_wire_mm = hit_to_wire_vector.Mag() / dd4hep::mm; // Explicitly cast to mm, no matter what the default unit is

            // Update the smallest distance to wire of all hits in the cell
            if (distance_to_wire_mm < closest_distance_to_wire_mm) {
                // Save the closest hit and related values
                closest_simhit = &simhit;
                closest_simhit_position = simhit_position;
                closest_hit_to_wire_vector = hit_to_wire_vector;
                closest_distance_to_wire_mm = distance_to_wire_mm;
            }
            
            // Integrate the energy deposited
            edep_sum += simhit.getEDep();

            /////////////////////////
            // CLUSTER INFORMATION //
            /////////////////////////

            auto mcparticle = simhit.getParticle();
            auto object_id = mcparticle.getObjectID();
            
            // Check if the particle is already in the map
            if (cluster_info_map.find(object_id) == cluster_info_map.end()) {
                // If not, create a new entry
                cluster_info_map[object_id] = ParticleClusterInfo();
                double mass = mcparticle.getMass();
                double momentum = edm4hep::utils::magnitude(mcparticle.getMomentum());
                // As simplfication: just take the betagamma value from the first hit we encounter
                double beta_gamma = momentum / mass;
                cluster_info_map[object_id].beta_gamma = beta_gamma;
                cluster_info_map[object_id].path_length = simhit.getPathLength();
            } else {
                // If the particle is already present, update the existing entry
                cluster_info_map[object_id].path_length += simhit.getPathLength();
            }

        }

        // Total number of clusters in cell built from each particle's contribution
        unsigned int total_nclusters = 0;

        // Calculate the number of clusters for each particle in the cell
        for (const auto& [object_id, particle_info] : cluster_info_map) {
            // Get number of clusters per length from delphes
            // Output from delphes function is in 1/m, so to convert to 1/mm we need to scale accordingly
            double nclusters_per_mm = m_delphesTrkUtil.Nclusters(particle_info.beta_gamma, m_GasSel.value()) / 1000.0;
            double nclusters_mean = nclusters_per_mm * particle_info.path_length;

            std::poisson_distribution<int> poisson_dist(nclusters_mean);
            total_nclusters += poisson_dist(random_engine);
        }



        ////////////////////////////////
        // POSITION AND TIME SMEARING //
        ////////////////////////////////

        double digihit_time;
        double digihit_distance_to_wire_mm;
        edm4hep::Vector3d digihit_position;

        if (closest_simhit) {
            // Do the digitisation based on the closest hit and the integrated edep and path length

            // TODO: update time with realisitc value
            digihit_time = closest_simhit->getTime();

            // xy smearing
            double smearing_xy = gauss_xy(random_engine);
            digihit_distance_to_wire_mm = std::max(0.0, closest_distance_to_wire_mm + smearing_xy);

            // z smearing
            double smearing_z = gauss_z(random_engine);
            auto hit_projection_on_the_wire = closest_simhit_position + closest_hit_to_wire_vector;
            TVector3 wire_direction_ez = (m_dch_info->Calculate_wire_vector_ez(layer, nphi)).Unit();
            hit_projection_on_the_wire += smearing_z * wire_direction_ez;

            // Convert to edm4hep vector and cast to mm
            digihit_position = this->Convert_TVector3_to_EDM4hepVector(hit_projection_on_the_wire, 1.0 / dd4hep::mm);



        } else {
            error() << "Could not find a closest hit for CellID: " << cellID << endmsg;
            continue;
        }

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

        sense_wire_hit.setCellID(cellID);
        sense_wire_hit.setType(0);
        sense_wire_hit.setQuality(0);
        sense_wire_hit.setTime(digihit_time);
        sense_wire_hit.setEDep(edep_sum);
        sense_wire_hit.setEDepError(0.0);
        sense_wire_hit.setPosition(digihit_position);
        sense_wire_hit.setPositionAlongWireError(m_z_resolution);
        sense_wire_hit.setWireAzimuthalAngle(WireAzimuthalAngle);
        sense_wire_hit.setWireStereoAngle(WireStereoAngle);
        sense_wire_hit.setDistanceToWire(digihit_distance_to_wire_mm);
        sense_wire_hit.setDistanceToWireError(m_xy_resolution);

        // Clusters are added to the SenseWireHit as vector member containing the number of electrons in each cluster
        // The length of this vector is the total number of clusters in the cell
        // For now, we do not calculate the size of each cluster, so we just fill it with a dummy value
        for (unsigned int i = 0; i < total_nclusters; ++i) {
            sense_wire_hit.addToNElectrons(-1);
        }

    }

    return output;

}


// DCHdigi_v02::get_drift_time(double distance_to_wire_mm) const {
//     // Calculate the drift time based on the distance to the wire
//     // This is a preliminary implementation that needs to be updated with a realistic model
//     // For now, we use a simple linear model with a constant drift velocity

//     double drift_velocity_um_per_ns = 30.0;

//     // Convert distance to wire from mm to um
//     double distance_to_wire_um = distance_to_wire_mm * 1000.0;

//     // Calculate the drift time in ns
//     double drift_time_ns = distance_to_wire_um / drift_velocity_um_per_ns;
//     return drift_time_ns;
// }


// DCHdigi_v02::get_signal_travel_time(double distance_to_readout_mm) const {
//     // Calculate the time it takes for the signal to travel along the wire to the readout electronics
//     // Assume 2/3 of the speed of light for the signal propagation speed

//     double speed_of_light_mm_per_ns = 299.792458;
//     double signal_speed_mm_per_ns = speed_of_light_mm_per_ns * (2.0 / 3.0);

//     // Calculate the signal travel time in ns
//     double signal_travel_time_ns = distance_to_readout_mm / signal_speed_mm_per_ns;
//     return signal_travel_time_ns;
// }
