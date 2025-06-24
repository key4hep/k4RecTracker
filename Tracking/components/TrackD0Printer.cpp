#include "Gaudi/Property.h"
#include "GaudiKernel/MsgStream.h"
#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "k4FWCore/Consumer.h"
#include "podio/RelationRange.h"
#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Type alias for TrackCollection to improve readability
using TrackColl = edm4hep::TrackCollection;
using TP = edm4hep::TrackParams;

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
    info() << std::string(25, '*') << endmsg;
    info() << "New event" << endmsg;
    info() << std::string(25, '*') << endmsg;

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

  float getSigmaPhi(const edm4hep::TrackState& ts) const { return std::sqrt(ts.getCovMatrix(TP::phi, TP::phi)); }

  std::string printValueUnc(const std::string& strTrType, const std::string& vName, const int maxLabelWidth,
                            const float varValue) const {
    return fmt::format("{:<10}{:>{}}{:>10.5f}", strTrType, vName, maxLabelWidth, varValue);
  }

  // Helper function to process each track (either SiTrack or CluTrack)
  void processTrack(const edm4hep::Track& track, const std::string& trackType) const {
    auto trackStates = track.getTrackStates(); // RelationRange<TrackState>

    // Check if there are enough track states (e.g., third track state)
    std::string varName = "phi";
    std::string sigmaVarName = "sigma " + varName;
    int maxVarWidth = std::max(varName.size(), sigmaVarName.size()) + 2;

    // Check if there are enough track states (e.g., third track state)
    if (trackStates.size() > 2) {
      const edm4hep::TrackState& state = trackStates[2]; // Get the third track state
      info() << printValueUnc(trackType, varName, maxVarWidth, state.phi) << endmsg;
      info() << printValueUnc(trackType, sigmaVarName, maxVarWidth, getSigmaPhi(state)) << endmsg;
    } else {
      warning() << trackType << " has less than 3 track states, skipping " + varName + " print" << endmsg;
    }
  }
};

// Declare the consumer component
DECLARE_COMPONENT(TrackD0Printer)