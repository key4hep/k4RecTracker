/*
 * MergeTrackerHitPlanes: merge several TrackerHitPlane collections (and their reco->sim
 * link collections) into one, so that one CombinatorialTrackFinder instance can find
 * tracks across several subdetectors that are digitised into separate collections
 * (e.g. vertex barrel + disks + silicon wrapper).
 *
 * The outputs are podio subset collections: they reference the original hits, so every
 * hit keeps its ObjectID (collectionID, index) and the reco->sim links stay valid.
 * Because the indices of the different input collections all start at 0, the track
 * finder keys hits on (collectionID, index).
 */
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include "k4FWCore/Transformer.h"

#include <string>
#include <tuple>
#include <vector>

struct MergeTrackerHitPlanes final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>(
              const std::vector<const edm4hep::TrackerHitPlaneCollection*>&,
              const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&)> {

  MergeTrackerHitPlanes(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {KeyValues{"InputHitCollections", {"VTXBDigis"}},
                          KeyValues{"InputRecoSimLinkCollections", {"VTXBSimDigiLinks"}}},
                         {KeyValues{"OutputHitCollection", {"MergedTrackerHits"}},
                          KeyValues{"OutputRecoSimLinkCollection", {"MergedTrackerHitSimLinks"}}}) {}

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& hitColls,
             const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& linkColls) const override {

    edm4hep::TrackerHitPlaneCollection hits;
    hits.setSubsetCollection(true);
    for (const auto* coll : hitColls) {
      if (!coll)
        continue;
      for (const auto hit : *coll)
        hits.push_back(hit);
    }

    edm4hep::TrackerHitSimTrackerHitLinkCollection links;
    links.setSubsetCollection(true);
    for (const auto* coll : linkColls) {
      if (!coll)
        continue;
      for (const auto link : *coll)
        links.push_back(link);
    }

    debug() << "Merged " << hits.size() << " hits from " << hitColls.size() << " collections and " << links.size()
            << " links from " << linkColls.size() << " collections" << endmsg;
    return std::make_tuple(std::move(hits), std::move(links));
  }
};

DECLARE_COMPONENT(MergeTrackerHitPlanes)
