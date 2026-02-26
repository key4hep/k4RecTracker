#include "utils.hpp"

#include "edm4hep/Track.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "podio/RelationRange.h"

#include "Gaudi/Property.h"
#include "GaudiKernel/MsgStream.h"

#include "k4FWCore/Consumer.h"

#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Type aliases for improved readability
using TrackColl = edm4hep::TrackCollection;
using TP = edm4hep::TrackParams;
using TS = edm4hep::TrackState;

// Consumer that processes track collections and prints phi values
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

    printInStars(this, "New Event", n_stars);

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
  Gaudi::Property<size_t> n_stars{this, "nStars", 20, "line-width of message in star box"};

  float getSigmaPhi(const edm4hep::TrackState& ts) const { return std::sqrt(ts.getCovMatrix(TP::phi, TP::phi)); }

  std::string printValueUnc(const std::string& strTrType, const std::string& vName, const int maxLabelWidth,
                            const float varValue) const {
    return fmt::format("{:<10}{:>{}}{:>10.5f}", strTrType, vName, maxLabelWidth, varValue);
  }

  // Helper function to process each track (either SiTrack or CluTrack)
  void processTrack(const edm4hep::Track& track, const std::string& trackType) const {

    std::string varName = "phi";
    std::string sigmaVarName = "sigma " + varName;
    int maxVarWidth = std::max(varName.size(), sigmaVarName.size()) + 2;

    // assuming there is only one Track State at IP
    if (auto trackAtIP = std::ranges::find(track.getTrackStates(), TS::AtIP, &TS::location);
        trackAtIP != track.getTrackStates().end()) {
      info() << printValueUnc(trackType, varName, maxVarWidth, trackAtIP->phi) << endmsg;
      info() << printValueUnc(trackType, sigmaVarName, maxVarWidth, getSigmaPhi(*trackAtIP)) << endmsg;
    } else {
      fatal() << fmt::format("No track at IP found for {}!", trackType) << endmsg;
    }
  }
};

// Declare the consumer component
DECLARE_COMPONENT(TrackD0Printer)
