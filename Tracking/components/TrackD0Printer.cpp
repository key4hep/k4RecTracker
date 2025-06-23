#include "Gaudi/Property.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "k4FWCore/Consumer.h"
#include "podio/RelationRange.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Type alias for TrackCollection to improve readability
using TrackColl = edm4hep::TrackCollection;

// Consumer that processes track collections and prints D0 values
struct TrackD0Printer final : k4FWCore::Consumer<void(const TrackColl&, const TrackColl&)> {
  // Constructor: define the input collections (ClupatraTracks and SiTracks)
  TrackD0Printer(const std::string& name, ISvcLocator* svcLoc)
      : Consumer(name, svcLoc,
                 {
                     KeyValues("InputSiTracks", {"SiTracks"}),
                     KeyValues("InputCluTracks", {"ClupatraTracks"}),
                 }) {}

  // This function will be called to process the data
  void operator()(const TrackColl& inSiTracks, const TrackColl& inCluTracks) const override {

    debug() << "Received SiTracks collection with " << inSiTracks.size() << " tracks" << endmsg;
    debug() << "Received ClupatraTracks collection with " << inCluTracks.size() << " tracks" << endmsg;

    // Check that both collections are of the same size
    if (inSiTracks.size() != inCluTracks.size()) {
      fatal() << "Track collections have different sizes: SiTracks (" << inSiTracks.size() << ") and ClupatraTracks ("
              << inCluTracks.size() << "). Exiting!" << endmsg;
      return;
    }

    // Process both collections (SiTracks and ClupatraTracks) at the same time
    for (size_t i = 0; i < inSiTracks.size(); ++i) {

      // Process SiTrack
      processTrack(inSiTracks[i], "SiTrack");

      // Process CluTrack
      processTrack(inCluTracks[i], "CluTrack");
    }
  }

private:
  // Optional: You can add a property to filter or adjust behavior if needed (e.g., track state index).
  Gaudi::Property<int> m_trackStateIndex{this, "TrackStateIndex", 2, "Index of track state to print (default 2)"};

  // Helper function to process each track (either SiTrack or CluTrack)
  void processTrack(const edm4hep::Track& track, const std::string& trackType) const {
    auto trackStates = track.getTrackStates(); // RelationRange<TrackState>

    // Check if there are enough track states (e.g., third track state)
    if (trackStates.size() > 2) {
      const edm4hep::TrackState& state = trackStates[2];          // Get the third track state
      std::cout << trackType << " D0: " << state.D0 << std::endl; // Print D0 value
    } else {
      warning() << trackType << " has less than 3 track states, skipping D0 print" << endmsg;
    }
  }
};

// Declare the consumer component
DECLARE_COMPONENT(TrackD0Printer)