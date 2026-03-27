/*
 * Copyright (c) 2014-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>

#include "Gaudi/Algorithm.h"

#include "k4FWCore/Transformer.h"

#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

/** @struct PerfectTrackFinder
 *
 *  @brief Algorithm that builds "perfect" tracks using MC truth information.
 *
 *  PerfectTrackFinder is a k4FWCore MultiTransformer that produces an
 *  edm4hep::TrackCollection starting from TrackerHit–SimTrackerHit link
 *  collections and the MCParticle collection.
 *
 *  For each MCParticle in the event, the algorithm:
 *
 *   - Creates one output edm4hep::Track.
 *   - Scans all provided TrackerHitSimTrackerHitLink collections
 *     (both planar and wire).
 *   - Selects the digitized tracker hits whose associated SimTrackerHit
 *     originates from the given MCParticle.
 *   - Adds those hits to the corresponding track.
 *
 *  The result is one "perfect" track per MCParticle, where hits are
 *  assigned using ground-truth Monte Carlo information instead of a
 *  pattern-recognition procedure.
 *
 *  This algorithm is typically used for:
 *   - Tracking performance studies
 *   - Validation of reconstruction algorithms
 *   - Efficiency and resolution benchmarking
 *
 *  Inputs:
 *   - InputPlanarHitCollections: vector of TrackerHitSimTrackerHitLink collections
 *   - InputWireHitCollections:   vector of TrackerHitSimTrackerHitLink collections
 *   - InputMCParticles:          MCParticle collection
 *
 *  Outputs:
 *   - OutputPerfectTracks: edm4hep::TrackCollection containing one MC-truth-based track per MCParticle
 *
 *  @author Andrea De Vita
 *  @date   2026-03
 *
 */

struct PerfectTrackFinder final : k4FWCore::MultiTransformer< std::tuple<edm4hep::TrackCollection> (
                                
                                const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,
                                const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,
                                const edm4hep::MCParticleCollection&)> {

  PerfectTrackFinder(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(name, svcLoc,
                         {

                            KeyValues("InputPlanarHitCollections", {"InputPlanarHitCollections"}),
                            KeyValues("InputWireHitCollections", {"InputWireHitCollections"}),
                            KeyValues("InputMCParticles", {"InputMCParticles"})

                         },
                         {

                            KeyValues("OutputPerfectTracks", {"OutputPerfectTracks"})

                         }) {}

  std::tuple<edm4hep::TrackCollection> 
    operator()(const std::vector<   const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& planarHitLinks,
                                    const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& wireHitLinks,
                                    const edm4hep::MCParticleCollection& mcParticles) const override {
    
                                        
    //////////////////////////////////
    ////////// PERFECT TRACKING //////
    //////////////////////////////////

    // Create a new TrackCollection for storing the output tracks
    edm4hep::TrackCollection outputTracks;

    // Loop over MCParticles to create perfect tracks
    for (const auto& mcParticle : mcParticles) {
    
        auto mcParticleObjectId = mcParticle.getObjectID();
        std::vector<std::pair<float, edm4hep::TrackerHit>> hitsWithTime;

        // Planar hits
        for (const auto& planarHitLinkCollection : planarHitLinks) {
            for (const auto& hitLink : *planarHitLinkCollection) {

                auto simHit = hitLink.getTo();
                auto digiHit = hitLink.getFrom();

                if (simHit.getParticle().getObjectID() == mcParticleObjectId) {
                    hitsWithTime.emplace_back(simHit.getTime(), digiHit);
                }
            }
        }

        // Wire hits
        for (const auto& wireHitLinkCollection : wireHitLinks) {
            for (const auto& hitLink : *wireHitLinkCollection) {

                auto simHit = hitLink.getTo();
                auto digiHit = hitLink.getFrom();

                if (simHit.getParticle().getObjectID() == mcParticleObjectId) {
                    hitsWithTime.emplace_back(simHit.getTime(), digiHit);
                }
            }
        }

        std::sort(
            hitsWithTime.begin(),
            hitsWithTime.end(),
            [](const auto& a, const auto& b) {
                return a.first < b.first;
            }
        );


        if (!hitsWithTime.empty()) {

            auto edm4hep_track = outputTracks.create();

            // Add all hits with their associated time
            for (const auto& [time, hit] : hitsWithTime) {
                edm4hep_track.addToTrackerHits(hit);
            }

            // Set track type as reconstructed
            edm4hep_track.setType(1);
        }

    }

    // Return the output collections as a tuple
    return std::make_tuple(std::move(outputTracks));
    
  }

};

DECLARE_COMPONENT(PerfectTrackFinder)