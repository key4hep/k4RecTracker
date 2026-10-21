#include "ElectronDecorator.h"

#include "utils.h"

#include <GaudiKernel/MsgStream.h>

#include <algorithm>

#include <tuple>

#include <cmath>

DECLARE_COMPONENT(ElectronDecorator)

// struct for storing cluster and assocated track from PFO
struct ClusterTrackMap {
    edm4hep::Cluster cluster;
    std::vector<edm4hep::Track> tracks;
};

// struct for storing brem Recovered cluster and assocated track for target PFO
struct bremRecoveredClusterTrackMap {
    std::vector<edm4hep::Cluster> clusters;
    std::vector<edm4hep::Track> tracks;
};

   

ElectronDecorator::ElectronDecorator( const std::string& name,ISvcLocator* svc)

    : Transformer(name,svc,{
              KeyValue{"InputRecoParticles", "PandoraPFOs"}
          },

          {
              KeyValue{"OutputRecoParticles",
                       "RecoveredPFOs"}
          })
{
std::cout << "pandoraReading constructor called" << std::endl;
}

edm4hep::ReconstructedParticleCollection 
ElectronDecorator::operator()(const edm4hep::ReconstructedParticleCollection& input) const {

    auto detector = m_geoSvc->getDetector();
    
    double eCalBarrelInnerR = detector->constantAsDouble("ECalBarrel_inner_radius")/dd4hep::mm;
    debug() << "eCalBarrelInnerR: " << eCalBarrelInnerR << endmsg;

    double eCalBarrelMaxZ = detector->constantAsDouble("ECalBarrel_half_length")/dd4hep::mm;
    debug() << "eCalBarrelMaxZ: " << eCalBarrelMaxZ << endmsg;

    double eCalEndCapInnerR = detector->constantAsDouble("ECalEndcap_inner_radius")/dd4hep::mm;
    debug() << "eCalEndCapInnerR: " << eCalEndCapInnerR << endmsg;

    double eCalEndCapOuterR = detector->constantAsDouble("ECalEndcap_outer_radius")/dd4hep::mm;
    debug() << "eCalEndCapOuterR: " << eCalEndCapOuterR << endmsg;

    double eCalEndCapInnerZ =detector->constantAsDouble("ECalEndcap_min_z")/dd4hep::mm;
    debug() << "eCalEndCapInnerZ: " <<eCalEndCapInnerZ << endmsg; 

    double m_Bz = getFieldFromCompact();

    // Forming the cluster->track map
    std::vector<ClusterTrackMap> clusterTrackMaps;

    std::vector<edm4hep::Cluster> tempClusters;

   //////////////////////////////////////////////////////////////////////////////////
   /////////////////Step1:loop through input PFOs ///////////////////////////////////
   /////////////////Step2:Unpack the PFOs to extract clusters and tracks ////////////
   /////////////////Step3:Make a map clusters: tracks ///////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////

    for (const auto& pfo : input) {  //PFO loop starts

        debug()
            << "PFO Energy = "
            << pfo.getEnergy()
            << " Tracks = "
            << pfo.getTracks().size()
            << " Clusters = "
            << pfo.getClusters().size()
            << endmsg;
        for (const auto& cluster : pfo.getClusters()) {  // PFO cluster loop starts

           if (!cluster.isAvailable())
              continue;
           
           tempClusters.push_back(cluster);           

           ClusterTrackMap entry;
           entry.cluster = cluster;

           for (const auto& track : pfo.getTracks()) { // PFO track loop starts
              if (track.isAvailable()) {
                  //auto state = track.getTrackStates(edm4hep::TrackState::AtIP);
                  edm4hep::TrackState  state;
                  auto stateHolder = track.getTrackState(edm4hep::TrackState::AtIP);

                  if (stateHolder.has_value()) {
                    state = stateHolder.value();
    
                  }

                  entry.tracks.push_back(track);
              }
           }  //  PFO track loop ends

           clusterTrackMaps.push_back(std::move(entry));
        }  //  PFO cluster loop starts ends

     } // PFO  loop ends

   //////////////////////////////////////////////////////////////////////////////////
   /////////////////Step 1:complete /////////////////////////////////////////////////
   /////////////////Step 2:complete /////////////////////////////////////////////////
   /////////////////Step 3:complete /////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////////////////
   /////////////////Step 4:Sort Cluster:Tracks map wrt cluster Energy ///////////////
   /////////////////Step 5:Sort Temporary Cluster Collection wrt cluster Energy /////
   //////////////////////////////////////////////////////////////////////////////////

     // sorting the cluster track map
     std::sort(clusterTrackMaps.begin(), clusterTrackMaps.end(),
      [](const ClusterTrackMap& a, const ClusterTrackMap& b) {
              return a.cluster.getEnergy() > b.cluster.getEnergy(); 
     });

     std::sort(tempClusters.begin(),tempClusters.end(),
       [](const edm4hep::Cluster& a, const edm4hep::Cluster& b){
            return a.getEnergy() > b.getEnergy();
       });

   //////////////////////////////////////////////////////////////////////////////////
   /////////////////Step 4:complete /////////////////////////////////////////////////
   /////////////////Step 5:complete /////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////////////////
   /////////////////Step 6:Extrapolate all tracks to clusters with utils function////
   /////////////////Step 7:Arrange the track wrt closest cluster to track distance //
   //////////////////////////////////////////////////////////////////////////////////

      for (auto& entry : clusterTrackMaps) {  //  Looping over clusterTrackMaps starts

          if (entry.tracks.size() <= 1)
            continue;

          const auto clusterPos = entry.cluster.getPosition();

          std::sort(entry.tracks.begin(), entry.tracks.end(),
            [&](const edm4hep::Track& a, const edm4hep::Track& b) {

               auto getDistance = [&](const edm4hep::Track& track) {

                 edm4hep::TrackState state=helixAtLastTrack(track, eCalBarrelInnerR, eCalBarrelMaxZ,  eCalEndCapInnerR, eCalEndCapOuterR,
                                             eCalEndCapInnerZ, m_Bz);

                 auto ref = state.referencePoint;
                 float dx = clusterPos.x - ref.x;
                 float dy = clusterPos.y - ref.y;
                 float dz = clusterPos.z - ref.z;

                 return std::sqrt(dx*dx + dy*dy + dz*dz);
               };

             return getDistance(a) < getDistance(b); // closest first
           });
        }  //Looping over clusterTrackMaps ends

   //////////////////////////////////////////////////////////////////////////////////
   /////////////////Step 6:complete /////////////////////////////////////////////////
   /////////////////Step 7:complete /////////////////////////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////


   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////Step 8: Application of brem interaction finding method (interaction of tangent from track curve to cluster)////
   ////////////8:a) Look for brem interaction point in EndCap Region  ///////////////////////////////////////////////////
   ////////////8:b) check brem interaction point in EndCap Region has x-y radius>barrel radius///////////////////////////
   ////////////8:c) if x-y radius>barrel radius look for brem interaction in barrel region //////////////////////////////
   ///////Step 9: Identify brem clusters from brem interaction point and its eta-phi margin /////////////////////////////
   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::map<int, std::vector<int>> clusterIndexMap;
        std::unordered_set<int> consumedClusters;
        std::vector<bremRecoveredClusterTrackMap> bremRecoveredClusterTrackMaps;
 

        for (const auto& entry : clusterTrackMaps) {  //  Looping over clusterTrackMaps starts for brem recovery

            const auto& cluster = entry.cluster;

            if (entry.tracks.empty()) {
               continue;
            }

            bool aStatus = false;
            const auto& track1st = entry.tracks[0];
            edm4hep::TrackState  trackState1st;
            auto state1st = track1st.getTrackState(edm4hep::TrackState::AtIP);

            if (state1st.has_value()) {
                trackState1st = state1st.value();
            }
            const auto refPt = ThreeVector(trackState1st.referencePoint.x, trackState1st.referencePoint.y, trackState1st.referencePoint.z);
            const auto momVec = getMomentum( trackState1st, m_Bz);
            double zPlane = eCalEndCapInnerZ;
            double radius = eCalBarrelInnerR;
            ThreeVector intersectZ;

            //!!! find brem iteraction in Endacap  region !!!//
            aStatus = GetPointInZ(zPlane,refPt,momVec,intersectZ);

            if (!aStatus) {
              throw std::runtime_error("BremRecoveryAlgorithm: failed to run BremRecoveryUtility::GetPointInZ!");
            }

            ThreeVector finalIntersectPt = intersectZ;

            if (radius <= 0.) { // this should never happen
              throw std::invalid_argument("BremRecoveryAlgorithm: ECAL barrel inner radius < 0!");
            } 

            //!!! if brem iteraction in Endacap  region has an x-y coordinate leading a radius greater  !!!//
            //!!! than barrel radius then look for brem interaction in barrel region !!!//
            
            if (intersectZ.magXY() > radius){ // the tangent intersects the barrel before reaching the endcap - recalculate it
              
              //!!! handled situation where beam interaction point at x(0),y(0),z(0) !!!//
              if(std::abs(refPt.x) < std::numeric_limits<float>::epsilon() &&
                 std::abs(refPt.y) < std::numeric_limits<float>::epsilon() &&
                 std::abs(refPt.z) < std::numeric_limits<float>::epsilon()){

                 aStatus = GetPointOnCircleForRefAtOrigin(radius,refPt,momVec,finalIntersectPt);

              }
              else{
                 aStatus = GetPointOnCircle(radius,refPt,momVec,finalIntersectPt);
              }
            }
            
           if (!aStatus) {
            throw std::runtime_error("BremRecoveryAlgorithm: failed to run BremRecoveryUtility::GetPointOnCircle!");
           }
           
           double dEtaMargin = m_dEtaMargin;
           double dPhiMargin = m_dPhiMargin;

           const auto& clusterPos1st = cluster.getPosition();
           const auto vecClus1 = ThreeVector(clusterPos1st.x,
                                            clusterPos1st.y,
                                            clusterPos1st.z);
           const double eta1st = vecClus1.eta();
           const double phi1st = vecClus1.phi();

           //!!! calculation of eta-phi margin for recoverd cluster with input parameters m_dEtaMargin and m_dPhiMargin !!!//
           const double maxEta = eta1st + dEtaMargin;
           const double minEta = eta1st - dEtaMargin;
           const double bremBandPhi = deltaPhi(finalIntersectPt.phi(),phi1st);
           const double maxDPhi = bremBandPhi > 0. ? bremBandPhi + dPhiMargin : dPhiMargin;
           const double minDPhi = bremBandPhi > 0. ? -dPhiMargin : bremBandPhi - dPhiMargin;

           //!!! looping over temporary clusters to check brem counter part of current cluster !!!//
           for (const auto& cluster2nd:tempClusters) { //  Looping over temporary clusters
              
             if (cluster2nd.getObjectID().index == cluster.getObjectID().index)  
               continue;
            
             if (!cluster2nd.isAvailable())
                continue;
             
             const auto& clusterPos2nd = cluster2nd.getPosition();
             const auto vecClus2 = ThreeVector(clusterPos2nd.x,
                                              clusterPos2nd.y,
                                              clusterPos2nd.z);
             const double eta2nd = vecClus2.eta();
             const double phi2nd = vecClus2.phi();
             const double finalDeltaPhi = deltaPhi(phi2nd,phi1st);

             if (minEta < eta2nd && eta2nd < maxEta && minDPhi < finalDeltaPhi && finalDeltaPhi < maxDPhi) {  //start of eta-phi margin inspection 
                debug()<<"===========found a brem candidate============"<<endmsg;
      
                int c1_idx = cluster.getObjectID().index;
                int c2_idx = cluster2nd.getObjectID().index;

                // Do not match if the secondary cluster has already been consumed as a key or value
                if (consumedClusters.find(c2_idx) == consumedClusters.end() && consumedClusters.find(c1_idx) == consumedClusters.end()) { //start if
          
                    clusterIndexMap[c1_idx].push_back(c2_idx);
        
                     // Mark the secondary cluster as consumed so it cannot be an entry elsewhere
                    consumedClusters.insert(c2_idx);
                }//end if
             }//end  of eta-phi margin inspection 
           } //  end Looping over temporary clusters
           int current_c1 = cluster.getObjectID().index;
           if (clusterIndexMap.find(current_c1) == clusterIndexMap.end() && consumedClusters.find(current_c1) == consumedClusters.end()) {
             clusterIndexMap[current_c1] = std::vector<int>(); // Empty value list
           }

        }  //  Looping over clusterTrackMaps ends for brem recovery

       //////////////////////////////////////////////////////////////////////////////////
       /////////////////Step 8:complete /////////////////////////////////////////////////
       /////////////////Step 9:complete /////////////////////////////////////////////////
       //////////////////////////////////////////////////////////////////////////////////

       //!!! Step 8 and 9 created a clusterIndexMap  of clusterID; key:electronClusterID and value: matched brem clusterID  !!!//
       //!!! using clusterIndexMap actual clusters and associated tracks are found from clusterTrackMaps in following steps !!!//

       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       //////////Step 10: Put brem recovered clusters associated tracks and remaining stuff in bremRecoveredClusterTrackMaps//////// 
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        if(clusterIndexMap.size()>0){
            std::set<int> recoveredTrackIndices;
            for (const auto& [primaryKey, matchedSecondaries] : clusterIndexMap) {
               
               bremRecoveredClusterTrackMap recoveredMap;
               recoveredTrackIndices.clear();
               
               for (const auto& ctm : clusterTrackMaps) {
         
                 int clusterID = ctm.cluster.getObjectID().index;
                 bool isPrimary = (clusterID == primaryKey);
                 bool isSecondary = (std::find(matchedSecondaries.begin(), matchedSecondaries.end(), clusterID) != matchedSecondaries.end()); 
                 if (isPrimary || isSecondary) {
                    recoveredMap.clusters.push_back(ctm.cluster);
                    for (const auto& track : ctm.tracks) {
                      int trackIndex = track.getObjectID().index;
                      if (recoveredTrackIndices.insert(trackIndex).second) {
                          recoveredMap.tracks.push_back(track);
                      }
                    }
                 }
               }
               bremRecoveredClusterTrackMaps.push_back(recoveredMap);
            }
        }
        for (const auto& ctm : clusterTrackMaps) {
           
            bremRecoveredClusterTrackMap remainingMap;
            bool isNotKey=true;
            bool isNotValue=true;
            int clusterID = ctm.cluster.getObjectID().index;
             
            if(clusterIndexMap.size()>0){
             
               isNotKey = (clusterIndexMap.find(clusterID) == clusterIndexMap.end());
               for (const auto& [key, matches] : clusterIndexMap) {
                   if (std::find(matches.begin(), matches.end(), clusterID) != matches.end()) {
                      isNotValue = false; 
                       break; // Found it as a value, so stop searching
                   }
               }
            }
            
            if (isNotKey && isNotValue) {
               remainingMap.clusters.push_back(ctm.cluster);
               for (const auto& track : ctm.tracks) {
                   remainingMap.tracks.push_back(track);

               }
               bremRecoveredClusterTrackMaps.push_back(std::move(remainingMap));
            }
        }

      ////////////////////////////////////////
      //////////Step 10: complete///////////// 
      ////////////////////////////////////////


      // building the output
      edm4hep::ReconstructedParticleCollection output;
      for (const auto& bremMap : bremRecoveredClusterTrackMaps) { // for bremMap starts

        // Create one ReconstructedParticle
        auto particle = output.create();
      
         float clusterEnergy = 0.0;

         // Add all associated tracks
         for (const auto& track : bremMap.tracks) { // for track starts 
            particle.addToTracks(track);
          }                                         // for track ends

         // Add all associated clusters
         for (const auto& cluster : bremMap.clusters) { // for cluster starts 
            particle.addToClusters(cluster);
            clusterEnergy += cluster.getEnergy();
         }                                                           // for cluster ends
          //std::cout << "clusterEnergy = " << clusterEnergy << std::endl;                                     
          particle.setEnergy(clusterEnergy);
       }// for bremMap ends
       return output;
}

