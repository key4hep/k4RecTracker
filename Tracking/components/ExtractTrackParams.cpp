#include "Gaudi/Property.h"
#include "edm4hep/TrackCollection.h"
#include "k4FWCore/Transformer.h"
#include "podio/UserDataCollection.h"
#include "printStars.h"

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
                             KeyValues("OutCollSiPhi", {"SiTrackPhi"}),
                             KeyValues("OutCollSiOmega", {"SiTrackOmega"}),
                             KeyValues("OutCollSiD0", {"SiTrackD0"}),
                             KeyValues("OutCollSiTanL", {"SiTrackTanL"}),
                             KeyValues("OutCollSiZ0", {"SiTrackZ0"}),
                             KeyValues("OutCollCluPhi", {"CluTrackPhi"}),
                             KeyValues("OutCollCluOmega", {"CluTrackOmega"}),
                             KeyValues("OutCollCluD0", {"CluTrackD0"}),
                             KeyValues("OutCollCluTanL", {"CluTrackTanL"}),
                             KeyValues("OutCollCluZ0", {"CluTrackZ0"}),

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
      siColls = processTrack(inSiTracks[i], "SiTrack");

      // Process CluTrack
      cluColls = processTrack(inCluTracks[i], "CluTrack");
    }

    return std::tuple_cat(std::move(siColls), std::move(cluColls));
  };

private:
  Gaudi::Property<size_t> n_stars{this, "nStars", 20, "line-width of message in star box"};

  std::tuple<FloatColl, FloatColl, FloatColl, FloatColl, FloatColl> processTrack(const edm4hep::Track& track,
                                                                                 const std::string& trackType) const {

    // assuming there is only one Track State at IP
    auto trackAtIP = std::ranges::find(track.getTrackStates(), TS::AtIP, &TS::location);
    if (trackAtIP != track.getTrackStates().end()) {
      verbose() << fmt::format("Track at IP found for {}.", trackType) << endmsg;
    } else {
      fatal() << fmt::format("No track at IP found for {}!", trackType) << endmsg;
    }
    FloatColl trackPhi;
    FloatColl trackOmega;
    FloatColl trackD0;
    FloatColl trackTanL;
    FloatColl trackZ0;
    trackPhi.push_back(trackAtIP->phi);
    trackOmega.push_back(trackAtIP->omega);
    trackD0.push_back(trackAtIP->D0);
    trackTanL.push_back(trackAtIP->tanLambda);
    trackZ0.push_back(trackAtIP->Z0);

    return std::make_tuple(std::move(trackPhi), std::move(trackOmega), std::move(trackD0), std::move(trackTanL),
                           std::move(trackZ0));
  }
};
DECLARE_COMPONENT(TrackParamExtractor)