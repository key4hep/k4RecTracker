#include "Gaudi/Property.h"
#include "edm4hep/TrackCollection.h"
#include "k4FWCore/Transformer.h"
#include "podio/UserDataCollection.h"
#include "printStars.h"
#include <edm4hep/Constants.h>

#include <optional>
#include <string>
#include <tuple>

// Which type of collection we are reading
using FloatColl = podio::UserDataCollection<float>;
using TrackColl = edm4hep::TrackCollection;
using TS = edm4hep::TrackState;

struct TrackParamExtractor final
    : k4FWCore::MultiTransformer<std::tuple<FloatColl, FloatColl, FloatColl, FloatColl, FloatColl, FloatColl, FloatColl,
                                            FloatColl, FloatColl, FloatColl>(const TrackColl&, const TrackColl&)> {
  TrackParamExtractor(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {
                             KeyValues("InputSiTracks", {"SiTracks"}),
                             KeyValues("InputCluTracks", {"ClupatraTracks"}),
                         },
                         {
                             KeyValues("OutCollSiD0", {"SiTrackD0"}),
                             KeyValues("OutCollSiPhi", {"SiTrackPhi"}),
                             KeyValues("OutCollSiOmega", {"SiTrackOmega"}),
                             KeyValues("OutCollSiZ0", {"SiTrackZ0"}),
                             KeyValues("OutCollSiTanL", {"SiTrackTanL"}),
                             KeyValues("OutCollCluD0", {"CluTrackD0"}),
                             KeyValues("OutCollCluPhi", {"CluTrackPhi"}),
                             KeyValues("OutCollCluOmega", {"CluTrackOmega"}),
                             KeyValues("OutCollCluZ0", {"CluTrackZ0"}),
                             KeyValues("OutCollCluTanL", {"CluTrackTanL"}),

                         }) {}

  // This is the function that will be called to transform the data
  // Note that the function has to be const, as well as the collections
  // we get from the input
  std::tuple<FloatColl, FloatColl, FloatColl, FloatColl, FloatColl, FloatColl, FloatColl, FloatColl, FloatColl,
             FloatColl>
  operator()(const TrackColl& inSiTracks, const TrackColl& inCluTracks) const override {

    printInStars(this, "New Event", n_stars);

    debug() << "Received SiTracks collection with " << inSiTracks.size() << " tracks" << endmsg;
    debug() << "Received ClupatraTracks collection with " << inCluTracks.size() << " tracks" << endmsg;

    // Check that both collections are of the same size
    if (inSiTracks.size() != inCluTracks.size()) {
      fatal() << "Track collections have different sizes: SiTracks (" << inSiTracks.size() << ") and ClupatraTracks ("
              << inCluTracks.size() << ")." << endmsg;
    }

    // Process both collections (SiTracks and ClupatraTracks) at the same time
    std::tuple<FloatColl, FloatColl, FloatColl, FloatColl, FloatColl> siColls;
    std::tuple<FloatColl, FloatColl, FloatColl, FloatColl, FloatColl> cluColls;
    for (size_t i = 0; i < inSiTracks.size(); ++i) {

      // Process SiTrack
      const auto oSiTrackStateIP = getOTrackAtIP(inSiTracks[i], "SiTrack");
      if (oSiTrackStateIP.has_value()) {
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::d0)>(siColls).push_back(oSiTrackStateIP->D0);
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::phi)>(siColls).push_back(oSiTrackStateIP->phi);
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::omega)>(siColls).push_back(oSiTrackStateIP->omega);
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::z0)>(siColls).push_back(oSiTrackStateIP->Z0);
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::tanLambda)>(siColls).push_back(
            oSiTrackStateIP->tanLambda);
      }

      // Process CluTrack
      const auto oCluTrackStateIP = getOTrackAtIP(inCluTracks[i], "CluTrack");
      if (oCluTrackStateIP.has_value()) {
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::d0)>(cluColls).push_back(oCluTrackStateIP->D0);
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::phi)>(cluColls).push_back(oCluTrackStateIP->phi);
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::omega)>(cluColls).push_back(oCluTrackStateIP->omega);
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::z0)>(cluColls).push_back(oCluTrackStateIP->Z0);
        std::get<static_cast<std::size_t>(edm4hep::TrackParams::tanLambda)>(cluColls).push_back(
            oCluTrackStateIP->tanLambda);
      }
    }

    return std::tuple_cat(std::move(siColls), std::move(cluColls));
  };

private:
  Gaudi::Property<size_t> n_stars{this, "nStars", 20, "line-width of message in star box"};

  std::optional<edm4hep::TrackState> getOTrackAtIP(const edm4hep::Track& track, const std::string& trackType) const {

    // assuming there is only one Track State at IP
    auto trackAtIP = std::ranges::find(track.getTrackStates(), TS::AtIP, &TS::location);
    if (trackAtIP != track.getTrackStates().end()) {
      verbose() << fmt::format("Track at IP found for {}.", trackType) << endmsg;
      return *trackAtIP;
    } else {
      fatal() << fmt::format("No track at IP found for {}!", trackType) << endmsg;
      return std::nullopt;
    }
  }
};
DECLARE_COMPONENT(TrackParamExtractor)