#include "DCHdigi_v02.h"

// edm4hep
#include "edm4hep/utils/vector_utils.h"

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
    m_event_counter = 0;
    return StatusCode::SUCCESS;
}

extension::SenseWireHitCollection DCHdigi_v02::operator()(const edm4hep::SimTrackerHitCollection& input,
                                                        const edm4hep::EventHeaderCollection& header) const {

    extension::SenseWireHitCollection output;
    
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
        const TVector3* closest_simhit_position = nullptr;
        const TVector3* closest_hit_to_wire_vector = nullptr;
        const double* closest_distance_to_wire = nullptr;

        double min_distance = std::numeric_limits<double>::max();
        double edep_sum = 0.0;
        double path_length_sum = 0.0;

        for (const auto& simhit : simhits) {
            // Need to find hit closest to the wire
            // Requires layer number, nphi, and hit position

            int layer = m_dch_info->CalculateILayerFromCellIDFields(
                m_decoder->get(cellID, "layer"),
                m_decoder->get(cellID, "superlayer")
            );
            int nphi = m_decoder->get(cellID, "nphi");
            // use dd4hep:mm as scale to convert into the dd4hep default unit system
            auto simhit_position = this->Convert_EDM4hepVector_to_TVector3(simhit.getPosition(), dd4hep::mm);

            auto hit_to_wire_vector = m_dch_info->Calculate_hitpos_to_wire_vector(layer, nphi, simhit_position);
            double distance_to_wire = hit_to_wire_vector.Mag();

            // Update the smallest distance to wire of all hits in the cell
            if (distance_to_wire < min_distance) {
                min_distance = distance_to_wire;
                closest_simhit = &simhit;
                closest_simhit_position = &simhit_position;
                closest_hit_to_wire_vector = &hit_to_wire_vector;
                closest_distance_to_wire = &distance_to_wire;
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
               << ", Min Time: " << min_distance 
               << ", EDep Sum: " << edep_sum 
               << ", Path Length Sum: " << path_length_sum
               << endmsg;

        double time;

        if (closest_simhit) {
            // Do the digitisation based on the closest hit and the integrated edep and path length
            time = closest_simhit->getTime();

            auto hit_projection_on_the_wire = (*closest_simhit_position) + (*closest_hit_to_wire_vector);


        } else {
            error() << "Could not find a closest hit for CellID: " << cellID << endmsg;
            continue;
        }

        sense_wire_hit.setCellID(cellID);
        sense_wire_hit.setType(0);
        sense_wire_hit.setQuality(0);
        sense_wire_hit.setTime(time);
        sense_wire_hit.setEDep(edep_sum);
        sense_wire_hit.setEDepError(path_length_sum);
        sense_wire_hit.setPosition(edm4hep::Vector3d(0.0, 0.0, 0.0));
        sense_wire_hit.setPositionAlongWireError(0.0);
        sense_wire_hit.setWireAzimuthalAngle(0.0);
        sense_wire_hit.setWireStereoAngle(0.0);
        sense_wire_hit.setDistanceToWire(min_distance);
        sense_wire_hit.setDistanceToWireError(0.0);
        


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