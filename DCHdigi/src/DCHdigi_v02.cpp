#include "DCHdigi_v02.h"

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

    // Create the output collection
    extension::SenseWireHitCollection output;
    
    debug() << "Processing event " << ++m_event_counter << endmsg;
    debug() << "Processing SimTrackerHitCollection with " << input.size() << " hits." << endmsg;
    std::unordered_map<uint64_t, std::vector<edm4hep::SimTrackerHit>> cellIDToHitsMap;
    // Loop over the inputs and save the cellIDs and corresponding hits to a map
    for (const auto& simhit: input){
        uint64_t cellID = simhit.getCellID();
        cellIDToHitsMap[cellID].push_back(simhit);
    }

    debug() << "Map contains " << cellIDToHitsMap.size() << " unique cellIDs." << endmsg;
    unsigned int counter = 0;
    for (const auto& [cellID, hits] : cellIDToHitsMap) {
        auto sense_wire_hit = output.create();
        sense_wire_hit.setCellID(cellID);
        sense_wire_hit.setType(0);
        sense_wire_hit.setQuality(0);

        float min_time = std::numeric_limits<float>::max();
        double edep_sum = 0.0;
        double path_length_sum = 0.0;
        for (const auto& hit : hits) {
            float time = hit.getTime();
            if (time < min_time) {
                min_time = time;
            }
            edep_sum += hit.getEDep();
            path_length_sum += hit.getPathLength();

            debug() << "Processing hit with CellID: " << cellID 
                   << ", Time: " << time 
                   << ", EDep: " << hit.getEDep() 
                   << ", PathLength: " << hit.getPathLength()
                   << endmsg;
        }

        debug() << "CellID: " << cellID 
               << ", Min Time: " << min_time 
               << ", EDep Sum: " << edep_sum 
               << ", Path Length Sum: " << path_length_sum
               << endmsg;
        sense_wire_hit.setTime(min_time);
        sense_wire_hit.setEDep(edep_sum);
        sense_wire_hit.setEDepError(path_length_sum);


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