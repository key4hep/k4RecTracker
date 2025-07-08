#include "DCHdigi_v02.h"

// Gaudi
#include <GaudiKernel/MsgStream.h>

// edm4hep
#include "edm4hep/utils/vector_utils.h"

// DD4hep
#include "DD4hep/Detector.h"

// STL
#include <random>

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

    unsigned int counter = 0;
    for (const auto& [cellID, simhits] : cell_map) {
        auto sense_wire_hit = output.create();

        // Save the closest hit to the wire
        // Digitised hit time and position will be based on closest hit
        const edm4hep::SimTrackerHit* closest_simhit = nullptr;
        // Save also the vectors that were calculated anyway to find the closest hit
        // They will be reused for digitisation of the hit
        TVector3 closest_simhit_position;
        TVector3 closest_hit_to_wire_vector;
        double closest_distance_to_wire;

        double min_distance = std::numeric_limits<double>::max();
        double edep_sum = 0.0;
        double path_length_sum = 0.0;

        // Finding hit closest to the wire requires layer number, nphi
        int layer = m_dch_info->CalculateILayerFromCellIDFields(
            m_decoder->get(cellID, "layer"),
            m_decoder->get(cellID, "superlayer")
        );
        int nphi = m_decoder->get(cellID, "nphi");

        for (const auto& simhit : simhits) {
            
            // Get hit position to calculate distance to wire
            // use dd4hep:mm as scale to convert into the dd4hep default unit system
            auto simhit_position = this->Convert_EDM4hepVector_to_TVector3(simhit.getPosition(), dd4hep::mm);

            auto hit_to_wire_vector = m_dch_info->Calculate_hitpos_to_wire_vector(layer, nphi, simhit_position);
            double distance_to_wire = hit_to_wire_vector.Mag();

            // Update the smallest distance to wire of all hits in the cell
            if (distance_to_wire < min_distance) {
                min_distance = distance_to_wire;

                // Save the closest hit and related values
                closest_simhit = &simhit;
                closest_simhit_position = simhit_position;
                closest_hit_to_wire_vector = hit_to_wire_vector;
                closest_distance_to_wire = distance_to_wire;
            }
            
            // Integrate the energy deposited and path length
            edep_sum += simhit.getEDep();
            path_length_sum += simhit.getPathLength();

            debug() << "Processing hit with CellID: " << cellID 
                    << ", Dist: " << distance_to_wire
                    << ", EDep: " << simhit.getEDep() 
                    << ", PathLength: " << simhit.getPathLength()
                    << endmsg;
        }

        debug() << "CellID: " << cellID 
               << ", Dist: " << min_distance 
               << ", EDep Sum: " << edep_sum 
               << ", Path Length Sum: " << path_length_sum
               << endmsg;

        double dighit_time;
        double dighit_distance_to_wire;
        edm4hep::Vector3d dighit_position;

        if (closest_simhit) {
            // Do the digitisation based on the closest hit and the integrated edep and path length

            // TODO: update time with realisitc value
            dighit_time = closest_simhit->getTime();

            // xy smearing
            double smearing_xy = gauss_xy(random_engine);
            dighit_distance_to_wire = std::max(0.0, closest_distance_to_wire + smearing_xy);

            // z smearing
            double smearing_z = gauss_z(random_engine);
            auto hit_projection_on_the_wire = closest_simhit_position + closest_hit_to_wire_vector;
            TVector3 wire_direction_ez = (this->m_dch_info->Calculate_wire_vector_ez(layer, nphi)).Unit();
            hit_projection_on_the_wire += smearing_z * wire_direction_ez;

            // Convert to edm4hep vector and cast to mm
            dighit_position = this->Convert_TVector3_to_EDM4hepVector(hit_projection_on_the_wire, 1.0 / dd4hep::mm);



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
/* END OF COPYING WIRE ANGLE CALCULATION FROM DCHdigi_v02 */

        sense_wire_hit.setCellID(cellID);
        sense_wire_hit.setType(0);
        sense_wire_hit.setQuality(0);
        sense_wire_hit.setTime(dighit_time);
        sense_wire_hit.setEDep(edep_sum);
        // For debugging, save the integrated path length as EDepError for now, to be updated for PR
        // IF YOU ARE DOING A CODE REVIEW AND EdepError IS STILL PATH LENGTH, PLEASE MAKE A COMMENT IN THE PR!
        sense_wire_hit.setEDepError(path_length_sum);
        sense_wire_hit.setPosition(dighit_position);
        sense_wire_hit.setPositionAlongWireError(m_z_resolution);
        sense_wire_hit.setWireAzimuthalAngle(WireAzimuthalAngle);
        sense_wire_hit.setWireStereoAngle(WireStereoAngle);
        sense_wire_hit.setDistanceToWire(min_distance);
        sense_wire_hit.setDistanceToWireError(m_xy_resolution);
        


        // Process each cellID and its corresponding hits
        // Here you would implement the logic to create SenseWireHits from the hits
        // For now, we just print the cellID and number of hits
        // debug() << "CellID: " << cellID << ", Number of Hits: " << hits.size() << endmsg;
        counter++;
    }
    debug() << "Processed " << counter << " unique cellIDs." << endmsg;


    /* Alvaro's code for reference 
    extension::MutableSenseWireHit oDCHdigihit;
    oDCHdigihit.setCellID(input_sim_hit.getCellID());
    oDCHdigihit.setType(type);
    oDCHdigihit.setQuality(quality);
    oDCHdigihit.setTime(input_sim_hit.getTime());
    oDCHdigihit.setEDep(input_sim_hit.getEDep());
    oDCHdigihit.setEDepError(eDepError);
    oDCHdigihit.setPosition(positionSW);
    oDCHdigihit.setPositionAlongWireError(m_z_resolution);
    oDCHdigihit.setWireAzimuthalAngle(WireAzimuthalAngle);
    oDCHdigihit.setWireStereoAngle(WireStereoAngle);
    oDCHdigihit.setDistanceToWire(distanceToWire);
    oDCHdigihit.setDistanceToWireError(m_xy_resolution);
    // For the sake of speed, let the dNdx calculation be optional
    if (m_calculate_dndx.value()) {
      auto [nCluster, nElectrons_v] = CalculateClusters(input_sim_hit, myRandom);
      // to return the total number of electrons within the step, do the following:
      //   int nElectronsTotal = std::accumulate( nElectrons_v.begin(), nElectrons_v.end(), 0);
      //   oDCHdigihit.setNElectronsTotal(nElectronsTotal);
      // to copy the vector of each cluster size to the EDM4hep data extension, do the following:
      for (auto ne : nElectrons_v)
        oDCHdigihit.addToNElectrons(ne);
    } */

    return output;

}
