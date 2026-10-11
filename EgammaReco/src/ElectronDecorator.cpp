#include "ElectronDecorator.h"

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

// struct for storing brem Recovered clustera and assocated track for target PFO
struct bremRecoveredClusterTrackMap {
    std::vector<edm4hep::Cluster> clusters;
    std::vector<edm4hep::Track> tracks;
};


struct ThreeVector {
    double x, y, z;
    
    ThreeVector() : x(0), y(0), z(0) {}
    ThreeVector(double x0, double y0, double z0) : x(x0), y(y0), z(z0) {}

    double mag() const { return std::sqrt(x*x + y*y + z*z); }
    double magXY() const { return std::sqrt(x*x + y*y); }
    double phi() const { return std::atan2(y,x); } // -pi to pi
    double theta() const { return std::acos(z / mag()); } // 0 to pi
    double eta() const { return -std::log(std::tan(theta()/2.)); }
  };

double deltaPhi(const double p1, const double p2) {
    double oo2pi = 1./(2.*M_PI);

    if (std::abs(p1-p2) <= M_PI)
      return p1 - p2;

    double n = std::round((p1-p2)*oo2pi);

    return p1 - p2 - n*2.*M_PI;
  }

ThreeVector getMomentum(const edm4hep::TrackState& state, double Bz)
{
    double phi        = state.phi;
    double omega      = state.omega;
    double tanLambda  = state.tanLambda;

    double pt = (0.299792458 * Bz / std::abs(omega))*dd4hep::MeV/dd4hep::GeV;

    double px = pt * std::cos(phi);
    double py = pt * std::sin(phi);
    double pz = pt * tanLambda;

    return ThreeVector(px, py, pz);
}

double getFieldFromCompact() {
    dd4hep::Detector& mainDetector = dd4hep::Detector::getInstance();
    const double position[3] = {0, 0, 0};      // position to calculate magnetic field at (the origin in this case)
    double magneticFieldVector[3] = {0, 0, 0}; // initialise object to hold magnetic field
    mainDetector.field().magneticField(position, magneticFieldVector); // get the magnetic field vector from DD4hep
    return magneticFieldVector[2] / dd4hep::tesla;                     // z component at (0,0,0)
  }


bool GetPointInZ(const double zPlane,        // z plane (endcap) coordinate
                   const ThreeVector& refPt,   // reference point
                   const ThreeVector& momVec,  // momentum
                   ThreeVector& intersectPt) { // output intersection point
    // protect from zero division
    if (std::abs(momVec.z) < std::numeric_limits<float>::epsilon()) {
      throw std::invalid_argument("BremRecoveryUtility: invalid momentum.z (zero division)");

      return false;
    }

    // check the assumption
    if (std::abs(refPt.z) > std::abs(zPlane)) {
      throw std::invalid_argument("BremRecoveryUtility: invalid refPoint.z (ref. point outside the endcap)");

      return false;
    }

    intersectPt.z = momVec.z > 0. ? std::abs(zPlane) : -std::abs(zPlane);
    //std::cout<<"intersectPt.z: "<<intersectPt.z<<std::endl;
    //std::cout<<"momVec.z: "<<momVec.z<<std::endl;
    const double genericTime = (intersectPt.z - refPt.z)/momVec.z;
    //std::cout<<"genericTime: "<<genericTime<<std::endl;
    //std::cout<<"momVec.x: "<<momVec.x<<std::endl;
    //std::cout<<"refPt.x: "<<refPt.x<<std::endl;
    intersectPt.x = refPt.x + momVec.x*genericTime;
    intersectPt.y = refPt.y + momVec.y*genericTime;

    return true;
  }

/* if track reherence point is not at origin */
bool GetPointOnCircle(const double radius,        // (barrel) radius
                        const ThreeVector& refPt,   // reference point
                        const ThreeVector& momVec,  // momentum
                        ThreeVector& intersectPt) { // output intersection point
    double refPtXY = refPt.magXY();
    double momVecXY = momVec.magXY();
    // protect from zero division
    if (std::abs(momVecXY) < std::numeric_limits<float>::epsilon()) {
      throw std::invalid_argument("BremRecoveryUtility: invalid pt (zero division)");

      return false;
    }

    // check the assumption
    if (std::abs(refPtXY) > radius) {
      throw std::invalid_argument("BremRecoveryUtility: invalid refPoint (ref. point outside the barrel)");

      return false;
    }
    // suppose a triangle OAP in the XY plane
    //      A        A is the reference point
    //      /\       O is the detector origin
    //     /  \      P is the intersection point
    //    ------     then the length OP = radius
    //   O      P

    // now do some middle school Euclidean geometry
    // cos(O+P) = refPt*momVec / |refPt|*|momVec|
    double normedInnerProd = (refPt.x*momVec.x + refPt.y*momVec.y) / (refPtXY*momVecXY);
    // sine rule : radius/sin(P+O) = |refPt|/sin(P)
    double sinP = refPtXY*std::sqrt(1. - normedInnerProd*normedInnerProd)/radius;
    double angOplusP = std::acos(normedInnerProd);
    double angP = std::asin(sinP); // always smaller than pi/2 by definition
    double angO = angOplusP - angP;
    // sine rule again : lenAP/sin(O) = |refPt|/sin(P)
    double lenAP = refPtXY*std::sin(angO)/sinP;
    double genericTime = lenAP/momVecXY; // this is what we wanted

    // now done with the geometry, come back to 3D
    intersectPt.x = refPt.x + momVec.x*genericTime;
    intersectPt.y = refPt.y + momVec.y*genericTime;
    intersectPt.z = refPt.z + momVec.z*genericTime;

    return true;
  }

/* if track reference at the origin */
bool GetPointOnCircleForRefAtOrigin(const double radius,        // (barrel) radius
                        const ThreeVector& refPt,   // reference point
                        const ThreeVector& momVec,  // momentum
                        ThreeVector& intersectPt) { // output intersection point
    double refPtXY = refPt.magXY();
    double momVecXY = momVec.magXY();
    // protect from zero division
    if (std::abs(momVecXY) < std::numeric_limits<float>::epsilon()) {
      throw std::invalid_argument("BremRecoveryUtility: invalid pt (zero division)");

      return false;
    }

    // check the assumption
    if (std::abs(refPtXY) > radius) {
      throw std::invalid_argument("BremRecoveryUtility: invalid refPoint (ref. point outside the barrel)");

      return false;
    }
    double genericTime = radius/momVecXY; // this is what we wanted

    // now done with the geometry, come back to 3D
    intersectPt.x = refPt.x + momVec.x*genericTime;
    intersectPt.y = refPt.y + momVec.y*genericTime;
    intersectPt.z = refPt.z + momVec.z*genericTime;

    return true;
  }

ElectronDecorator::ElectronDecorator( const std::string& name,ISvcLocator* svc)

    : Transformer(name,svc,{
              KeyValue{"InputRecoParticles", "PandoraPFOs"},
              KeyValue{"InputClusters", "PandoraClusters"},
              KeyValue{"InputTracks", "SiTracks_Refitted"}
          },

          {
              KeyValue{"OutputRecoParticles",
                       "RecoveredPFOs"}
          })
{
std::cout << "pandoraReading constructor called" << std::endl;
}

edm4hep::ReconstructedParticleCollection 
ElectronDresser::operator()(const edm4hep::ReconstructedParticleCollection& input,
                           const edm4hep::ClusterCollection& clusters,
                           const edm4hep::TrackCollection& tracks) const {

    auto detector = m_geoSvc->getDetector();
    
    double eCalBarrelInnerR = detector->constantAsDouble("ECalBarrel_inner_radius")/dd4hep::mm;
     //std::cout << "eCalBarrelInnerR: " << eCalBarrelInnerR << std::endl;

    double eCalBarrelMaxZ = detector->constantAsDouble("ECalBarrel_half_length")/dd4hep::mm;
     //std::cout << "eCalBarrelMaxZ: " << eCalBarrelMaxZ << std::endl;

    double eCalEndCapInnerR = detector->constantAsDouble("ECalEndcap_inner_radius")/dd4hep::mm;
     //std::cout << "eCalEndCapInnerR: " << eCalEndCapInnerR << std::endl;

    double eCalEndCapOuterR = detector->constantAsDouble("ECalEndcap_outer_radius")/dd4hep::mm;
     //std::cout << "eCalEndCapOuterR: " << eCalEndCapOuterR << std::endl;

    double eCalEndCapInnerZ =detector->constantAsDouble("ECalEndcap_min_z")/dd4hep::mm;
     //std::cout << "eCalEndCapInnerZ: " <<eCalEndCapInnerZ << std::endl; 

    double m_Bz = getFieldFromCompact();

    // Forming the cluster->track map
    std::vector<ClusterTrackMap> clusterTrackMaps;

    std::vector<edm4hep::Cluster> tempClusters;

    for (const auto& pfo : input) {  //PFO loop starts

        info()
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
                  auto state = track.getTrackStates(edm4hep::TrackState::AtIP);
                  const auto momVec = getMomentum(state, m_Bz);

                  float p = std::sqrt(momVec.x*momVec.x + momVec.y*momVec.y + momVec.z*momVec.z);

                  // massless approximation

                  Total_Input_PFO_Track_Energy = Total_Input_PFO_Track_Energy+p;
                  entry.tracks.push_back(track);
              }
           }  //  PFO track loop ends

           clusterTrackMaps.push_back(std::move(entry));
        }  //  PFO cluster loop starts ends

     } // PFO  loop ends

     // sorting the cluster track map
     std::sort(clusterTrackMaps.begin(), clusterTrackMaps.end(),
      [](const ClusterTrackMap& a, const ClusterTrackMap& b) {
              return a.cluster.getEnergy() > b.cluster.getEnergy(); 
     });

     std::sort(tempClusters.begin(),tempClusters.end(),
       [](const edm4hep::Cluster& a, const edm4hep::Cluster& b){
            return a.getEnergy() > b.getEnergy();
       });

      for (auto& entry : clusterTrackMaps) {  //  Looping over clusterTrackMaps starts

          if (entry.tracks.size() <= 1)
            continue;

          std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!there are more than one track==============="<<std::endl;

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

        std::vector<int> clusterIndex;
        std::vector<bremRecoveredClusterTrackMap> bremRecoveredClusterTrackMaps;

        for (const auto& entry : clusterTrackMaps) {  //  Looping over clusterTrackMaps starts for brem recovery

            const auto& cluster = entry.cluster;

            //std::cout << "Cluster energy = " << cluster.getEnergy() << std::endl;

            if (entry.tracks.empty()) {
             //std::cout << "  No associated tracks." << std::endl;
             continue;
            }

            if (!entry.tracks[0].isAvailable())
                continue;

            bool aStatus = false;
            const auto& track1st = entry.tracks[0];

            const auto& trackState1st = track1st.getTrackStates(edm4hep::TrackState::AtIP);
            const auto refPt = ThreeVector(trackState1st.referencePoint.x, trackState1st.referencePoint.y, trackState1st.referencePoint.z);
            const auto momVec = getMomentum( trackState1st, m_Bz);
            double zPlane = eCalEndCapInnerZ;
            double radius = eCalBarrelInnerR;
            ThreeVector intersectZ;
            aStatus = GetPointInZ(zPlane,refPt,momVec,intersectZ);

            if (!aStatus) {
              throw std::runtime_error("BremRecoveryAlgorithm: failed to run BremRecoveryUtility::GetPointInZ!");
            }

            ThreeVector finalIntersectPt = intersectZ;

            if (radius <= 0.) { // this should never happen
              throw std::invalid_argument("BremRecoveryAlgorithm: ECAL barrel inner radius < 0!");
            } 
            
            if (intersectZ.magXY() > radius){ // the tangent intersects the barrel before reaching the endcap - recalculate it
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
           const double maxEta = eta1st + dEtaMargin;
           const double minEta = eta1st - dEtaMargin;
           const double bremBandPhi = deltaPhi(finalIntersectPt.phi(),phi1st);
           const double maxDPhi = bremBandPhi > 0. ? bremBandPhi + dPhiMargin : dPhiMargin;
           const double minDPhi = bremBandPhi > 0. ? -dPhiMargin : bremBandPhi - dPhiMargin;

           std::vector<edm4hep::Cluster> bremCands;
           for (const auto& cluster2nd:tempClusters) { //  Looping over temporary clusters
              
             if (cluster2nd.getObjectID().index == cluster.getObjectID().index){  // start if 
               continue;
             }      // end if 
             if (!cluster2nd.isAvailable())
                continue;
             
             const auto& clusterPos2nd = cluster2nd.getPosition();
             const auto vecClus2 = ThreeVector(clusterPos2nd.x,
                                              clusterPos2nd.y,
                                              clusterPos2nd.z);
             const double eta2nd = vecClus2.eta();
             const double phi2nd = vecClus2.phi();
             const double finalDeltaPhi = deltaPhi(phi2nd,phi1st);

             if (minEta < eta2nd && eta2nd < maxEta && minDPhi < finalDeltaPhi && finalDeltaPhi < maxDPhi) {  //start if
                std::cout<<"===========found a brem candidate============"<<std::endl;
      
                if (std::find(clusterIndex.begin(), clusterIndex.end(), cluster2nd.getObjectID().index) == clusterIndex.end()){ //start if
                    bremCands.push_back(cluster2nd);
                    clusterIndex.push_back(cluster2nd.getObjectID().index);
                }//end if
        
                if (std::find(clusterIndex.begin(), clusterIndex.end(), cluster.getObjectID().index) == clusterIndex.end()){ //start if
                    bremCands.push_back(cluster);
                    clusterIndex.push_back(cluster.getObjectID().index);
        
                }//end if
             }//end if
           } //  end Looping over temporary clusters
        }  //  Looping over clusterTrackMaps ends for brem recovery

        // for all selected and recovered clusters

        bremRecoveredClusterTrackMap recoveredMap;

        // Keep track IDs here so that tracks are not duplicated
        std::unordered_set<int> recoveredTrackIndices;

        for (const auto& ctm : clusterTrackMaps) {
           
            int clusterID = ctm.cluster.getObjectID().index;

           // Is this cluster one of the selected clusters?
            if (std::find(clusterIndex.begin(), clusterIndex.end(), clusterID)!= clusterIndex.end()){

             // Store the cluster
              recoveredMap.clusters.push_back(ctm.cluster);

             // Store associated tracks, but only once
              for (const auto& track : ctm.tracks) {
                //std::cout << "=========== Inside CTM tracks brem ============ " << std::endl;
                int trackIndex = track.getObjectID().index;
                //std::cout << "ctm.tracks = " << trackIndex << std::endl;

                if (recoveredTrackIndices.insert(trackIndex).second) {
                   recoveredMap.tracks.push_back(track);
                }
              }
           }
        }

        // --------------------------------------------------
       //  Add the recovered/selected group
       // --------------------------------------------------

       if (!recoveredMap.clusters.empty()) {
        bremRecoveredClusterTrackMaps.push_back(
        std::move(recoveredMap)
        );
       }

      // --------------------------------------------------
      //       Store every remaining cluster separately
      // --------------------------------------------------

       for (const auto& ctm : clusterTrackMaps) {

         int clusterID = ctm.cluster.getObjectID().index;

         // Skip clusters already included in recoveredMap
         if (std::find(clusterIndex.begin(), clusterIndex.end(), clusterID)!= clusterIndex.end()) {
           continue;
         }

         bremRecoveredClusterTrackMap remainingMap;

         // One cluster
         remainingMap.clusters.push_back(ctm.cluster);

         // Its associated tracks
         for (const auto& track : ctm.tracks) {
          //std::cout << "=========== Inside CTM tracks============ " << std::endl;
          //std::cout <<"track.getObjectID().index:::"<<track.getObjectID().index<<std::endl;
          remainingMap.tracks.push_back(track);
         }

         bremRecoveredClusterTrackMaps.push_back(
          std::move(remainingMap)
        );
      }

      // building the output
      edm4hep::ReconstructedParticleCollection output;
      for (const auto& bremMap : bremRecoveredClusterTrackMaps) { // for bremMap starts

        // Create one ReconstructedParticle
        auto particle = output.create();
      
         float clusterEnergy = 0.0;
         float pfoEnergy = 0.0;

         // Add all associated tracks
         for (const auto& track : bremMap.tracks) { // for track starts 
            particle.addToTracks(track);
          }                                         // for track ends

         // Add all associated clusters
         for (const auto& cluster : bremMap.clusters) { // for cluster starts 
            particle.addToClusters(cluster);
            clusterEnergy += cluster.getEnergy();
         }                                                           // for cluster ends
          std::cout << "clusterEnergy = " << clusterEnergy << std::endl;                                     
          particle.setEnergy(clusterEnergy);
       }// for bremMap ends


}
