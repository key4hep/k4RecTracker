/*
 * Copyright (c) 2020-2024 Key4hep-Project.
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

#include "GenfitTrack.h"
#include "KalmanFitterInfo.h"
#include "AbsFitterInfo.h"
#include "EventDisplay.h"


namespace GenfitInterface {

    GenfitTrack::GenfitTrack(   const edm4hep::Track& track,
                                const bool skipTrackOrdering,
                                const dd4hep::rec::DCH_info* dch_info, const dd4hep::DDSegmentation::BitFieldCoder* decoder)

        :   
            m_originalTrack(track),
            m_posInit(0., 0., 0.), 
            m_momInit(0., 0., 0.), 
            m_covInit(6),
            m_genfitTrackRep(nullptr), 
            m_genfitTrack(nullptr), 
            m_edm4hepTrack(),
        
            m_dch_info(dch_info),
            m_dc_decoder(decoder)

            
    {   

        CheckInitialization();
        OrderHits(track, skipTrackOrdering);

    }

    GenfitTrack::~GenfitTrack() {}

    /**
    * @brief Check if required Genfit components are properly initialized.
    * 
    * This method verifies whether the singleton instances of `genfit::FieldManager`
    * and `genfit::MaterialEffects` have been initialized. These components are essential
    * for tracking operations in the Genfit framework. 
    * 
    * If either component is not initialized, an error message is printed to `std::cerr` 
    * and the program exits with `EXIT_FAILURE`.
    * 
    * @note This method should be called before performing any tracking-related operations 
    * that depend on magnetic field or material effects.
    */
    void GenfitTrack::CheckInitialization() {

        if (!genfit::FieldManager::getInstance()->isInitialized()) {
            std::cerr << "Error: FieldManager is not initialized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (!genfit::MaterialEffects::getInstance()->isInitialized()) {
            std::cerr << "Error: MaterialEffects is not initialized!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

    }

    /**
    * @brief Orders the tracker hits of a track according to their spatial progression.
    *
    * This method determines a suitable starting point among the track hits,
    * typically the one closest to the interaction region (|z| minimum).
    * In case of nearly transverse tracks (almost perpendicular to the z axis),
    * a fallback criterion based on spatial radius is applied.
    *
    * All hits are then sorted according to their 3D distance from the chosen
    * starting point, ensuring a consistent ordering along the track trajectory.
    * The ordered hits are stored in the internal MutableTrack object for
    * subsequent processing (e.g. track fitting).
    *
    * @param track The input track whose tracker hits must be reordered.
    * @param skipTrackOrdering If true, skip the track ordering step.
    */
    void GenfitTrack::OrderHits(const edm4hep::Track& track, bool skipTrackOrdering) {   

        if (skipTrackOrdering) {
            
            for (const auto& hit : track.getTrackerHits()) {
                
                m_edm4hepTrack.addToTrackerHits(hit);
            }
            return;
        }

        const auto hits = track.getTrackerHits();
        std::vector<std::pair<float, std::size_t>> distIndex;

        // Track endpoints in cylindrical radius
        double rMin = std::numeric_limits<double>::max();
        double firstX = 0., firstY = 0., firstZ = 0.;

        // Loop over hits to find the hit with smallest r = sqrt(x^2 + y^2)
        for (const auto& hit : hits) {
            const auto p = hit.getPosition();
            double r = std::hypot(p.x, p.y);

            if (r < rMin) {
                rMin = r;
                firstX = p.x;
                firstY = p.y;
                firstZ = p.z;
            }
        }


        // Compute distances from the chosen starting point
        distIndex.reserve(hits.size());

        for (std::size_t i = 0; i < hits.size(); ++i) {
            const auto p = hits[i].getPosition();
            const float d = std::hypot(p.x - firstX, p.y - firstY, p.z - firstZ);

            distIndex.emplace_back(d, i);
        }

        // Sort hits along the track
        std::ranges::sort(distIndex, {}, &std::pair<float, std::size_t>::first);

        m_edm4hepTrack = edm4hep::MutableTrack();
        for (const auto& [_, idx] : distIndex) {
            m_edm4hepTrack.addToTrackerHits(hits[idx]);

        }
        
    }

    /**
    * @brief Initialize the GenfitTrack object with configurable initialization strategies.
    *
    * This method initializes the internal state of the `GenfitTrack` object,
    * setting the initial position, momentum direction, and covariance matrix
    * according to the selected initialization strategy.
    *
    * The initialization requires the magnetic field component Bz and supports
    * optional hit limiting and custom initialization parameters.
    *
    * The following steps are performed:
    *
    * 1. Optional hit reduction:
    *    - If LimitHits is true, the number of tracker hits is reduced using
    *      the `LimitNumberHits(Epsilon, Window)` method.
    *
    * 2. Track parameter initialization (controlled by InitializationType):
    *
    *    - InitializationType == 0 (default):
    *        - Uses the first two tracker hits.
    *        - m_posInit is set to the position of the first hit.
    *        - m_momInit is set to the unit vector from the first to the second hit.
    *        - The covariance matrix is computed via `ComputeInitialCovarianceMatrix`.
    *
    *    - InitializationType == 1 (refined):
    *        - Calls `ComputeInitialParameters(Bz)` to obtain improved
    *          estimates of the initial position and momentum.
    *        - The covariance matrix is computed accordingly.
    *
    *    - InitializationType == 2 (track state):
    *        - Extracts initial parameters from a specified track state location.
    *        - The position is taken from the track state's reference point.
    *        - The momentum is computed from the track state's helix parameters
    *          using the magnetic field Bz and physical constants.
    *
    *    - InitializationType == 3 (custom):
    *        - Uses user-provided Init_position and Init_momentum.
    *        - Both parameters must be specified.
    *        - The covariance matrix is computed accordingly.
    *
    * @param Bz Magnetic field component along the z-axis.
    * @param LimitHits If true, reduces the number of hits used for initialization.
    * @param InitializationType Defines the initialization strategy:
    *        0 = default (first two hits),
    *        1 = refined (computed parameters),
    *        2 = track state (parameters from a specific track state location),
    *        3 = custom (user-defined position and momentum).
    * @param TrackStateLocation Required if InitializationType == 2; specifies the track state location to extract parameters from.
    * @param Init_position Optional custom initial position (required if InitializationType == 3).
    * @param Init_momentum Optional custom initial momentum (required if InitializationType == 3).
    * @param Epsilon Optional parameter for hit limiting (required if LimitHits == true).
    * @param Window Optional parameter for hit limiting (required if LimitHits == true).
    *
    * @note The method terminates execution if required parameters are missing
    *       or if an unknown InitializationType is specified.
    */
    void GenfitTrack::InitializeTrack(double Bz, bool LimitHits, int InitializationType, std::optional<int> TrackStateLocation, std::optional<TVector3> Init_position, std::optional<TVector3> Init_momentum, std::optional<double> Epsilon, std::optional<int> Window) {

        double c_mm_s = 2.998e11;
        double a = 1e-15 * c_mm_s;

        if (LimitHits)
        {   
            if (Epsilon.has_value() && Window.has_value())  
            {
                LimitNumberHits(Epsilon.value(), Window.value()); 
            }
            else
            {
                std::cerr << "Missing Values for Epsilon or Window!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }

        HelperInitialization InitInfo = ComputeInitialParameters(Bz);
        m_charge_hypothesis = InitInfo.Charge;

        // Inizialization
        if (InitializationType == 0) // default
        {

            auto hits_for_genfit = m_edm4hepTrack.getTrackerHits();
            TVector3 firstHit_referencePoint;
            TVector3 secondHit_referencePoint;
            int index_loopHit = 0;
            for (auto hit : hits_for_genfit) {

                if (index_loopHit == 0) {
        
                    auto position = hit.getPosition();
                    firstHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z); // cm

                }
                if (index_loopHit == 1) {

                    auto position = hit.getPosition();
                    secondHit_referencePoint = TVector3(dd4hep::mm * position.x, dd4hep::mm * position.y, dd4hep::mm * position.z); // cm
                
                }

                index_loopHit++;
                if (index_loopHit > 1) break;

            }

            m_posInit = firstHit_referencePoint;                                       // cm
            m_momInit = (secondHit_referencePoint - firstHit_referencePoint).Unit();   // GeV/c
            m_covInit = ComputeInitialCovarianceMatrix(Bz, a);
            
        }
        else if (InitializationType == 1) // refined 
        {

            m_posInit = InitInfo.Position;
            m_momInit = InitInfo.Momentum;
            m_covInit = ComputeInitialCovarianceMatrix(Bz, a);


        }
        else if (InitializationType == 2) // take from trackState
        {

            if (TrackStateLocation.has_value())
            {
            

                auto trackState = m_originalTrack.getTrackStates();
                for (auto ts : trackState) {

                    if (ts.location == TrackStateLocation.value())
                    {   
                        m_VP_referencePoint = TVector3(ts.referencePoint.x * dd4hep::mm, ts.referencePoint.y * dd4hep::mm, ts.referencePoint.z * dd4hep::mm ); // cm
                        
                        double pT = a * std::abs(Bz) / std::abs(ts.omega);  // GeV/c
                        double px = pT * std::cos(ts.phi);
                        double py = pT * std::sin(ts.phi);
                        double pz = pT * ts.tanLambda;

                        m_posInit = m_VP_referencePoint; // cm
                        m_momInit = TVector3(px, py, pz); // GeV/c

                        auto covMatrixHelix = ts.covMatrix;
                        std::array<double,25> covMatrixHelix_array{};
                        for (int i = 0; i < 5; ++i) {
                            for (int j = 0; j < 5; ++j) {
                                covMatrixHelix_array[i*5 + j] = static_cast<double>(covMatrixHelix[i*6 + j]);
                            }
                        }
                        auto covMatrixCartesian_array = CovarianceMatrixHelixToCartesian(covMatrixHelix_array, m_posInit, m_momInit, m_VP_referencePoint, Bz);
                        TMatrixDSym covMatrixCartesian(6);
                        for (int i = 0; i < 6; ++i) {
                            for (int j = 0; j <= i; ++j) {
                                covMatrixCartesian(i, j) = covMatrixCartesian_array[i*6 + j];
                                if (i != j) {
                                    covMatrixCartesian(j, i) = covMatrixCartesian_array[i*6 + j]; // Symmetric
                                }
                            }
                        }
                        m_covInit = covMatrixCartesian; 

                        break;
                    }

                }
            }
            else
            {
                std::cerr << "Missing Value for TrackStateLocation!" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        else if (InitializationType == 3) // custom 
        {

            if (Init_position.has_value() && Init_momentum.has_value())  
            {
                m_posInit = Init_position.value();
                m_momInit = Init_momentum.value();
                m_covInit = ComputeInitialCovarianceMatrix(Bz, a);
            }
            else
            {
                std::cerr << "Missing Values for Init_position or Init_momentum!" << std::endl;
                std::exit(EXIT_FAILURE);
            }

        }
        else
        {   
            std::cerr << "Unknown InitializationType!" << std::endl;
            std::exit(EXIT_FAILURE);
        }

    }

    /**
    * @brief Limits the number of tracker hits in the track based on a local extremum 
    *        along the y-coordinate after smoothing.
    *
    * This method identifies the first significant local extremum (maximum or minimum) 
    * in the y-coordinate of the track hits, after applying a moving average smoothing. 
    * It then truncates the track to include only hits up to that extremum.
    *
    * Procedure:
    * 1. Extracts the z and y positions of all tracker hits.
    * 2. Applies a moving average smoothing to the y-coordinates using the specified window.
    * 3. Computes the first derivative of the smoothed y with respect to z.
    * 4. Finds the first index where the derivative changes sign and exceeds the 
    *    specified threshold epsilon, indicating a local extremum.
    * 5. If such an extremum is found, rebuilds the track keeping only hits up to that index.
    *
    * @param epsilon       Threshold for detecting significant changes in the derivative.
    * @param smoothWindow  Size of the smoothing window (number of hits) for the moving average.
    */
    void GenfitTrack::LimitNumberHits(double epsilon, int smoothWindow) {

        auto hits = m_edm4hepTrack.getTrackerHits();
        int maxHit = hits.size();
        int n = maxHit;

        if (maxHit == 0)
        {
            std::cerr << "Internal edm4hep::Track is empty." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (n < smoothWindow || n < 3) maxHit = -1;
        int W = smoothWindow / 2;

        std::vector<double> z_val;
        std::vector<double> y_val_raw;
        z_val.reserve(hits.size());
        y_val_raw.reserve(hits.size());

        // Extract positions
        for (const auto& hit : hits) {
            auto pos = hit.getPosition();
            z_val.push_back(pos.z);    
            y_val_raw.push_back(pos.y);
        }

        // Apply moving average smoothing
        std::vector<double> y_val;
        y_val.reserve(hits.size());
        for (std::size_t i = 0; i < hits.size(); ++i) {
            int start = (i < static_cast<std::size_t>(W)) ? 0 : i - W;
            int end   = std::min(i + W + 1, hits.size());

            double sum = 0.0;
            for (int j = start; j < end; ++j) sum += y_val_raw[j];
            y_val.push_back(sum / (end - start));
        }

        // Compute first derivative
        std::vector<double> d(n, 0.0);
        d[0] = (y_val[1] - y_val[0]) / (z_val[1] - z_val[0]);
        for (std::size_t i = 1; i < static_cast<std::size_t>(W - 1); ++i) {
            double dz = z_val[i+1] - z_val[i-1];
            d[i] = (dz != 0.0) ? (y_val[i+1] - y_val[i-1]) / dz : 0.0;
        }
        d[n-1] = (y_val[n-1] - y_val[n-2]) / (z_val[n-1] - z_val[n-2]);

        // Find first significant extremum
        bool MaxFound = false;
        for (int i = 1; i < n - 1; ++i) {
            if ((d[i-1] >  epsilon && d[i] < -epsilon) ||  
                (d[i-1] < -epsilon && d[i] >  epsilon))    
            {
                maxHit = i;
                MaxFound = true;
                break;
            }
        }

        if (!MaxFound) maxHit = -1;

        // Rebuild track keeping only hits up to maxHit
        if (maxHit != -1)
        {
            auto temp_track = m_edm4hepTrack;
            m_edm4hepTrack = edm4hep::MutableTrack();

            int idx_fill_track = 0;
            auto temp_hits = temp_track.getTrackerHits();
            for (const auto& hit : temp_hits) {
                m_edm4hepTrack.addToTrackerHits(hit);
                ++idx_fill_track;
                if (idx_fill_track >= maxHit) break;
            }
        }
    }

    /**
    * @brief Computes the initial 6x6 covariance matrix in Cartesian coordinates.
    *
    * The method performs the following steps:
    *
    * - Initializes a 5x5 covariance matrix in helix parameter space
    *   (d0, phi, omega, z0, tanLambda), with diagonal elements:
    *     - d0          : (0.05 cm)^2
    *     - phi         : (0.1 rad)^2
    *     - omega       : (0.5·omega)^2
    *     - z0          : (0.1·z0)^2
    *     - tanLambda   : (0.1)^2
    * - Computes the curvature parameter:
    *       omega = |a * Bz / pt|
    *   where pt is derived from m_momInit.
    * - Converts the covariance matrix from helix to Cartesian
    *   representation (6x6) using CovarianceMatrixHelixToCartesian().
    * - Fills a TMatrixDSym (6x6) with the converted values.
    * - Returns the Cartesian state covariance matrix.
    *
    * @param Bz Magnetic field component along the z-axis.
    * @param a  Conversion factor used in curvature calculation.
    *
    * @return 6x6 symmetric covariance matrix in Cartesian coordinates.
    */
    TMatrixDSym GenfitTrack::ComputeInitialCovarianceMatrix(double Bz, double a) {
        std::array<double,25> C_helix{};
        C_helix.fill(0.0);

        // columns: 0=d0, 1=phi, 2=omega, 3=z0, 4=tanLambda

        C_helix[0*5 + 0] = 0.05 * 0.05;                         // d0 : 500 um = 0.05 cm
        C_helix[1*5 + 1] = 0.1 * 0.1;                           // phi : 0.1 rad
             
        double pt = m_momInit.Perp();
        double omega = std::abs(a * Bz / pt * dd4hep::mm);      // a * Bz / pt

        C_helix[2*5 + 2] = std::pow(0.5 * omega, 2);            // omega : 0.5*omega 
        C_helix[3*5 + 3] = std::pow(0.1 * m_posInit.Z(), 2);    // z0 : 0.1*z0 
        C_helix[4*5 + 4] = 0.1 * 0.1;                           // tanLambda : 0.1  

        std::array<double,36> C_cart = 
            CovarianceMatrixHelixToCartesian(
                C_helix,
                m_posInit,
                m_momInit,
                m_VP_referencePoint,
                Bz
            );
        
        TMatrixDSym covState(6);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                covState(i,j) = C_cart[i*6 + j];
            }
        }

        return covState;

    }


    GenfitTrack::HelperInitialization GenfitTrack::ComputeInitialParameters(double Bz) {

        struct Point3D
        {
            double x;
            double y;
            double z;
        };

        // FIT CIRCLE TO XY PROJECTION
        std::vector<Point3D> points;
        auto hits = m_edm4hepTrack.getTrackerHits();
        for (const auto& hit : hits) {
            
            auto p = hit.getPosition();

            points.push_back(Point3D(p.x, p.y, p.z));
        }

        std::vector<Point2D_xy> points_xy;
        for (const auto& p : points) {
            points_xy.push_back(Point2D_xy(p.x, p.y));
        }

        FastCircleFit circle(points_xy);
        Point2D_xy closestPoint = circle.closestPointTo(points_xy[0]);
        Point2D_xy tangent_xy = circle.tangentAtPCA(closestPoint, points_xy[1]);

        double rho = circle.rho();
        double init_pT = std::abs(rho * 0.3 * Bz) / 1000; 

        TVector3 init_mom = TVector3(tangent_xy.x * init_pT, tangent_xy.y * init_pT, 0);
        
        // FIT Z
        double pZ = 0;

        size_t N = points.size();

        std::vector<Point2D_xy> points_Rz;
        for (const auto& p : points) {

            double R = std::sqrt(std::pow(p.x,2) + std::pow(p.y,2));
            double z_coord = p.z;
            points_Rz.push_back(Point2D_xy(R, z_coord));
        }

        double sumR = 0.0;
        double sumZ = 0.0;
        double sumRZ = 0.0;
        double sumR2 = 0.0;

        for (const auto& p : points_Rz) {
            sumR  += p.x;
            sumZ  += p.y;
            sumRZ += p.x * p.y;
            sumR2 += p.x * p.x;
        }

        double numerator   = N * sumRZ - sumR * sumZ;
        double denominator = N * sumR2 - sumR * sumR;

        pZ = numerator / denominator * init_pT;
        init_mom.SetZ(pZ);

        double a = numerator / denominator;
        sumZ = 0.0;
        for (const auto& p : points_Rz) sumZ += p.y;
        double b = (sumZ - a * sumR) / N;

        double z_PCA = b;

        // CHARGEs
        double x0_ = circle.x0();
        double y0_ = circle.y0();
        TVector3 B_field = TVector3(0., 0., 2.);

        TVector3 pos(closestPoint.x, closestPoint.y, z_PCA);
        TVector3 center(x0_, y0_, 0);

        TVector3 curvDir = center - pos;
        TVector3 lorentzDir = init_mom.Cross(B_field);
        lorentzDir.SetZ(0);
        curvDir.SetZ(0);

        int charge = (lorentzDir.Dot(curvDir) > 0) ? +1 : -1;

        // Fill helpers
        HelperInitialization helper;
        helper.Position = TVector3(closestPoint.x * dd4hep::mm, closestPoint.y * dd4hep::mm, z_PCA * dd4hep::mm);
        helper.Momentum = init_mom;
        helper.Charge = charge;

        return helper;

    };


    /**
    * @brief Builds the genfit::Track from the initialized state and tracker hits.
    *
    * The method performs the following steps:
    *
    * - Sets the signed particle hypothesis (including charge sign handling
    *   for leptons).
    * - Deletes any existing genfit::TrackRep and genfit::Track objects.
    * - Constructs the 6D state vector (x, y, z, px, py, pz) from m_posInit
    *   and m_momInit.
    * - Creates a genfit::RKTrackRep using the signed particle hypothesis.
    * - Instantiates a genfit::Track with the track representation,
    *   state vector, and m_covInit.
    * - Iterates over tracker hits stored in m_edm4hepTrack.
    * - Converts each hit into a Genfit measurement:
    *     - edm4hep::TrackerHitPlane: planar measurement (detID = 0)
    *     - edm4hep::SenseWireHit: wire measurement   (detID = 1)
    * - Wraps each measurement into a genfit::TrackPoint
    *   and inserts it into the genfit::Track.
    * - Terminates execution if an unsupported hit type is found.
    *
    * @param particle_hypotesis PDG code hypothesis (sign adjusted internally).
    * @param debug_lvl Debug verbosity level passed to measurement constructors.
    *
    * @note InitializeTrack() must be called before this method.
    */
    void GenfitTrack::CreateGenFitTrack(int particle_hypotesis, int debug_lvl) {

        m_signed_particle_hypothesis = particle_hypotesis;
        if (particle_hypotesis == 11 || particle_hypotesis == 13)
        {
            m_signed_particle_hypothesis = - m_charge_hypothesis * particle_hypotesis;

        }
        else
        {   
            m_signed_particle_hypothesis = m_charge_hypothesis * particle_hypotesis;

        }   

        delete m_genfitTrackRep;
        delete m_genfitTrack;

        // Create stateVec
        TVectorD stateVec(6);
        
        stateVec[0] = m_posInit.X();
        stateVec[1] = m_posInit.Y();
        stateVec[2] = m_posInit.Z();

        stateVec[3] = m_momInit.X();
        stateVec[4] = m_momInit.Y();
        stateVec[5] = m_momInit.Z();


        m_genfitTrackRep = new genfit::RKTrackRep(m_signed_particle_hypothesis);
        m_genfitTrack = new genfit::Track(m_genfitTrackRep, stateVec, m_covInit);

        auto hits_for_genfit = m_edm4hepTrack.getTrackerHits();
 
        int hit_idx(0);
        int detID(-1);
        for (auto hit : hits_for_genfit)
        {

            auto cellID0 = hit.getCellID();
            if (hit.isA<edm4hep::TrackerHitPlane>())
            {
                
                detID = 0;
                auto planar_hit =  hit.as<edm4hep::TrackerHitPlane>();
                GenfitInterface::PlanarMeasurement measurement = GenfitInterface::PlanarMeasurement(planar_hit, detID, ++hit_idx, debug_lvl);
                m_genfitTrack->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack));
                
            }
            else if (hit.isA<edm4hep::SenseWireHit>()) 
            {
                
                detID = 1;
                auto wire_hit =  hit.as<edm4hep::SenseWireHit>();
                GenfitInterface::WireMeasurement measurement = GenfitInterface::WireMeasurement(wire_hit, m_dch_info, m_dc_decoder, detID, ++hit_idx, debug_lvl);
                m_genfitTrack->insertPoint(new genfit::TrackPoint(measurement.getGenFit(), m_genfitTrack));

                       
            } 
            else 
            {

                std::cerr << "Error: No hits with cellID: " << cellID0 << std::endl;
                std::exit(EXIT_FAILURE);
             
            }

        }

    }

    /**
    * @brief Fit the Genfit track using the Deterministic Annealing Filter (DAF).
    * 
    * This function performs a forward and backward track fit using the Genfit DAF algorithm,
    * and populates the EDM4hep track states at the IP, first hit, and last hit positions.
    * 
    * @param Beta_init  Initial annealing parameter (temperature).
    * @param Beta_final Final annealing parameter.
    * @param Beta_steps Number of steps in the annealing schedule.
    * @param Bz         Value of z-component of the magnetic field at the center of the detector
    * @param debug_lvl  debug level: output if > 0
    * 
    * @return true if the fit was successful, false otherwise.
    */
    bool GenfitTrack::Fit(double Beta_init = 100., double Beta_final=0.1, double Beta_steps=10, double Bz = 2.0, int debug_lvl = 0) {

        edm4hep::Track Track_temp = m_edm4hepTrack;
        for (size_t i = 0; i < Track_temp.trackStates_size(); ++i) {

            Track_temp.getTrackStates(i) = edm4hep::TrackState();
        }

        try{

            // Initialize the genfit fitter
            genfit::DAF* m_genfitFitter = new genfit::DAF(true, 1e-3,1e-3);

            m_genfitFitter->setAnnealingScheme(Beta_init,Beta_final,Beta_steps);
            m_genfitFitter->setProbCut(1e-5);
            m_genfitFitter->setConvergenceDeltaWeight(1e-2); 

            int debug_lvl_fit = 0;
            if (debug_lvl > 2) debug_lvl_fit = 1;
            m_genfitFitter->setDebugLvl(debug_lvl_fit);

            // Process track
            genfit::Track genfitTrack = *m_genfitTrack;
            genfit::AbsTrackRep* trackRep = genfitTrack.getTrackRep(0);
            m_genfitFitter->processTrackWithRep(&genfitTrack, trackRep);

            // Update edm4hep track state
            genfit::MeasuredStateOnPlane fittedState;
            TVector3 gen_position, gen_momentum;
            TMatrixDSym covariancePosMom(6);
            
            double x_ref;
            double y_ref;
            double z_ref; 

            double pz;  
            double pt;

            double d0;
            double z0;
            double phi;
            double omega;  
            double tanLambda; 

            double c_mm_s = 2.998e11;
            double a = 1e-15 * c_mm_s;
            
            if (m_genfitFitter->isTrackFitted(&genfitTrack, trackRep))
            {

            
                // trackState First Hit
                fittedState = genfitTrack.getFittedState();
                fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                auto stateVecFirstHit = fittedState.getState();
                
                edm4hep::TrackState trackStateFirstHit;
                x_ref = gen_position.X();     // cm
                y_ref = gen_position.Y();     // cm
                z_ref = gen_position.Z();     // cm
                pz = gen_momentum.Z();        // gev
                pt = gen_momentum.Perp();     // gev
                
                auto infoComputeD0Z0_firstHit = PCAInfo(gen_position, gen_momentum, m_charge_hypothesis, m_VP_referencePoint, Bz);

                d0 = (( - (m_VP_referencePoint.X() - infoComputeD0Z0_firstHit.PCA.X()) ) * sin(infoComputeD0Z0_firstHit.Phi0) + 
                        (m_VP_referencePoint.Y() - infoComputeD0Z0_firstHit.PCA.Y())*cos(infoComputeD0Z0_firstHit.Phi0) ) / dd4hep::mm;   // mm
                z0 = ( infoComputeD0Z0_firstHit.PCA.Z() - m_VP_referencePoint.Z() ) / dd4hep::mm;                                // mm
                phi = gen_momentum.Phi();                                                                               // rad    

                tanLambda = pz / pt;
                omega =  std::abs(a * Bz / pt);
                if (m_charge_hypothesis < 0) omega = - omega;

                trackStateFirstHit.D0 = d0;
                trackStateFirstHit.Z0 = z0;
                trackStateFirstHit.phi = phi;
                trackStateFirstHit.omega = omega;
                trackStateFirstHit.tanLambda = tanLambda;
                trackStateFirstHit.time = 0.;

                trackStateFirstHit.referencePoint = edm4hep::Vector3f(x_ref / dd4hep::mm, y_ref / dd4hep::mm, z_ref / dd4hep::mm);
                trackStateFirstHit.location = edm4hep::TrackState::AtFirstHit;

                // trackState lastHit
                fittedState = genfitTrack.getFittedState(genfitTrack.getNumPoints()-1);
                fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                auto stateVecLastHit = fittedState.getState();
                
                edm4hep::TrackState trackStateLastHit;
                x_ref = gen_position.X();       // cm
                y_ref = gen_position.Y();       // cm
                z_ref = gen_position.Z();       // cm
                pz = gen_momentum.Z();          // gev
                pt = gen_momentum.Perp();       // gev

                auto infoComputeD0Z0_lastHit = PCAInfo(gen_position, gen_momentum, m_charge_hypothesis, m_VP_referencePoint, Bz);

                d0 = (( - (m_VP_referencePoint.X() - infoComputeD0Z0_lastHit.PCA.X()) ) * sin(infoComputeD0Z0_lastHit.Phi0) + 
                        (m_VP_referencePoint.Y() - infoComputeD0Z0_lastHit.PCA.Y())*cos(infoComputeD0Z0_lastHit.Phi0) ) / dd4hep::mm;   // mm
                z0 = ( infoComputeD0Z0_lastHit.PCA.Z() - m_VP_referencePoint.Z() ) / dd4hep::mm;                                // mm
                phi = gen_momentum.Phi();                                                                               // rad    

                tanLambda = pz / pt;
                omega =  std::abs(a * Bz / pt);
                if (m_charge_hypothesis < 0) omega = - omega;

                trackStateLastHit.D0 = d0;
                trackStateLastHit.Z0 = z0;
                trackStateLastHit.phi = phi;
                trackStateLastHit.omega = omega;
                trackStateLastHit.tanLambda = tanLambda;
                trackStateLastHit.time = 0.;

                trackStateLastHit.referencePoint = edm4hep::Vector3f(x_ref / dd4hep::mm, y_ref / dd4hep::mm, z_ref / dd4hep::mm);
                trackStateLastHit.location = edm4hep::TrackState::AtLastHit;

                //take first fitted point
                genfit::TrackPoint* tp = genfitTrack.getPointWithFitterInfo(0);
                if (tp == NULL) {std::cout << "Track has no TrackPoint with fitterInfo! (but fitstatus ok?)"<<std::endl;}
                auto* fi = static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(trackRep));

                // extrapolate rep to target plane (IP)
                try{

                    fittedState=fi->getFittedState(true);
                    trackRep->extrapolateToLine(fittedState,TVector3(0,0,0),TVector3(0,0,1));

                    fittedState.getPosMomCov(gen_position, gen_momentum, covariancePosMom);
                    auto stateVecIP = fittedState.getState();

                    gen_momentum.SetX(-gen_momentum.X());
                    gen_momentum.SetY(-gen_momentum.Y());
                    gen_momentum.SetZ(-gen_momentum.Z());
                    
                    edm4hep::TrackState trackStateIP;
                    x_ref = gen_position.X();       // cm
                    y_ref = gen_position.Y();       // cm
                    z_ref = gen_position.Z();       // cm
                    pz = gen_momentum.Z();          // gev
                    pt = gen_momentum.Perp();       // gev

                    auto infoComputeD0Z0_IP = PCAInfo(gen_position, gen_momentum, m_charge_hypothesis, m_VP_referencePoint, Bz);

                    d0 = (( - (m_VP_referencePoint.X() - infoComputeD0Z0_IP.PCA.X()) ) * sin(infoComputeD0Z0_IP.Phi0) + 
                            (m_VP_referencePoint.Y() - infoComputeD0Z0_IP.PCA.Y())*cos(infoComputeD0Z0_IP.Phi0) ) / dd4hep::mm;     // mm
                    z0 = ( infoComputeD0Z0_IP.PCA.Z() - m_VP_referencePoint.Z() ) / dd4hep::mm;                                     // mm
                    phi = gen_momentum.Phi();                                                                                       // rad    

                    tanLambda = pz / pt;
                    omega =  std::abs(a * Bz / pt);
                    if (m_charge_hypothesis < 0) omega = - omega;

                    trackStateIP.D0 = d0;
                    trackStateIP.Z0 = z0;
                    trackStateIP.phi = phi;
                    trackStateIP.omega = omega;
                    trackStateIP.tanLambda = tanLambda;
                    trackStateIP.time = 0.;

                    trackStateIP.referencePoint = edm4hep::Vector3f(x_ref / dd4hep::mm, y_ref / dd4hep::mm, z_ref / dd4hep::mm);
                    trackStateIP.location = edm4hep::TrackState::AtIP;

                    if (debug_lvl > 0)
                    {

                        std::cout << "GenfitTrackFitter    DEBUG : TrackState at IP: " << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  D0: " << trackStateIP.D0 << " mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  Z0: " << trackStateIP.Z0 << " mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  phi: " << trackStateIP.phi << " rad" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  omega: " << trackStateIP.omega << " 1/mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  tanLambda: " << trackStateIP.tanLambda << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  location: " << trackStateIP.location << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  reference point: (" << trackStateIP.referencePoint.x << ", " << trackStateIP.referencePoint.y << ", " << trackStateIP.referencePoint.z << ") mm" << std::endl;
                        
                        std::cout << "\nGenfitTrackFitter    DEBUG : TrackState at First Hit: " << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  D0: " << trackStateFirstHit.D0 << " mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  Z0: " << trackStateFirstHit.Z0 << " mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  phi: " << trackStateFirstHit.phi << " rad" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  omega: " << trackStateFirstHit.omega << " 1/mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  tanLambda: " << trackStateFirstHit.tanLambda << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  location: " << trackStateFirstHit.location << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  reference point: (" << trackStateFirstHit.referencePoint.x << ", " << trackStateFirstHit.referencePoint.y << ", " << trackStateFirstHit.referencePoint.z << ") mm" << std::endl;

                        std::cout << "\nGenfitTrackFitter    DEBUG : TrackState at Last Hit: " << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  D0: " << trackStateLastHit.D0 << " mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  Z0: " << trackStateLastHit.Z0 << " mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  phi: " << trackStateLastHit.phi << " rad" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  omega: " << trackStateLastHit.omega << " 1/mm" << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  tanLambda: " << trackStateLastHit.tanLambda << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  location: " << trackStateLastHit.location << std::endl;
                        std::cout << "GenfitTrackFitter    DEBUG :  reference point: (" << trackStateLastHit.referencePoint.x << ", " << trackStateLastHit.referencePoint.y << ", " << trackStateLastHit.referencePoint.z << ") mm\n" << std::endl;

                    }   

                    m_edm4hepTrack.addToTrackStates(trackStateIP);
                    m_edm4hepTrack.addToTrackStates(trackStateFirstHit);
                    m_edm4hepTrack.addToTrackStates(trackStateLastHit);
                    
                }catch(...){

                    return false;
                }

                if (m_genfitFitter->isTrackFitted(&genfitTrack, trackRep))
                {
                    
                    m_edm4hepTrack.setChi2(genfitTrack.getFitStatus()->getChi2());
                    m_edm4hepTrack.setNdf(genfitTrack.getFitStatus()->getNdf());

                }
                else
                {

                    m_edm4hepTrack.setChi2(-1);
                    m_edm4hepTrack.setNdf(-1);
                    return false;

                }

            }
            else
            {

                m_edm4hepTrack.setChi2(-1);
                m_edm4hepTrack.setNdf(-1);
                return false;

            }

            return m_genfitFitter->isTrackFitted(&genfitTrack, trackRep);
            
        }
        catch(...)
        {
            return false;
        }
 
    }

    /**
    * @brief Propagation of the covariance matrix from helix-parameter representation to Cartesian coordinate representation
    *
    * @param C_helix        Covariance Matrix in helix-basis
    * @param Position_cm    Inizial Position in cm
    * @param Momentum_gev   Initial Momentum in gev/c
    * @param RefPoint_cm    Reference position (e.g. IP)
    * @param Bz             Bz
    * 
    * @return Covariance matrix in cartesian-basis
    */
    std::array<double, 36> GenfitTrack::CovarianceMatrixHelixToCartesian(   const std::array<double,25>& C_helix, TVector3 Position_cm, TVector3 Momentum_gev,  TVector3 RefPoint_cm, double Bz) {                                                       
                                                                
        double c_mm_s = 2.998e11;
        double a = 1e-15 * c_mm_s;

        // Assumption: the initialization refers to the PCA 
        double x_PCA = Position_cm.X();
        double y_PCA = Position_cm.Y();

        double px = Momentum_gev.X();
        double py = Momentum_gev.Y();
        double pz = Momentum_gev.Z();  

        double pt   = Momentum_gev.Perp(); 
        double phi0 = std::atan2(py, px);               

        double d0 = - (RefPoint_cm.X() - x_PCA) * sin(phi0) + (RefPoint_cm.Y() - y_PCA) * cos(phi0); // cm                                            // rad
        double tanLambda = pz / pt;                                                                  
        double omega = (std::abs(a * Bz / pt)) * dd4hep::mm;                                         // cm-1

        // Jacobian J[row][col] = d(cartesian variable) / d(helix parameter)
        // rows:    0=x, 1=y, 2=z, 3=px, 4=py, 5=pz
        // columns: 0=d0, 1=phi, 2=omega, 3=z0, 4=tanLambda

        double J[6][5] = {0};

        ///////
        // x //
        ///////

        J[0][0] = 1.0 / sin(phi0);  
        // dx / dd0

        J[0][1] = cos(phi0) / (sin(phi0) * sin(phi0)) * d0 - (y_PCA - RefPoint_cm.Y()) / (cos(phi0) * cos(phi0) * std::tan(phi0) * std::tan(phi0));
        // dx / dphi0

        J[0][2] =  0.0;  
        // dx / domega

        J[0][3] =  0.0;  
        // dx / dz0

        J[0][4] =  0.0;  
        // dx / dtanLambda


        ///////
        // y //
        ///////

        J[1][0] = - 1.0 / cos(phi0);  
        // dy / dd0

        J[1][1] =  - sin(phi0) / (cos(phi0) * cos(phi0)) * d0 - (RefPoint_cm.X() - x_PCA) / (cos(phi0) * cos(phi0));
        // dy / dphi0

        J[1][2] =  0.0;  
        // dy / domega

        J[1][3] =  0.0;  
        // dy / dz0

        J[1][4] =  0.0;  
        // dy / dtanLambda


        ///////
        // z //
        ///////
        
        J[2][0] =  0.0;  
        // dz / dd0

        J[2][1] =  0.0;  
        // dz / dphi0

        J[2][2] =  0.0;  
        // dz / domega

        J[2][3] =  1.0;  
        // dz / dz0

        J[2][4] =  0.0;  
        // dz / dtanLambda


        /////////////////////////////////////////////////////////////////////////////////////
        // px = p * cos(phi0) * sin(theta) = pT * cos(phi0) = a * |Bz| / omega * cos(phi0) //
        /////////////////////////////////////////////////////////////////////////////////////

        J[3][0] =  0.0;  
        // dpx / dd0

        J[3][1] = - px * std::tan(phi0);
        // dpx / dphi0

        J[3][2] = - px / std::abs(omega);
        // dpx / domega

        J[3][3] =  0.0;  
        // dpx / dz0

        J[3][4] = - Momentum_gev.Mag() * std::cos(phi0) * tanLambda / std::pow(1.0 + tanLambda * tanLambda, 1.5);
        // dpx / dtanLambda


        /////////////////////////////////////////////////////////////////////////////////////
        // py = p * sin(phi0) * sin(theta) = pT * sin(phi0) = a * |Bz| / omega * sin(phi0) //
        /////////////////////////////////////////////////////////////////////////////////////

        J[4][0] =  0.0;  
        // dpy / dd0

        J[4][1] = - py / std::tan(phi0);
        // dpy / dphi0

        J[4][2] = - py / std::abs(omega);
        // dpy / domega

        J[4][3] =  0.0;  
        // dpy / dz0

        J[4][4] = - Momentum_gev.Mag() * std::sin(phi0) * tanLambda / std::pow(1.0 + tanLambda * tanLambda, 1.5);
        // dpy / dtanLambda

        ////////////////////////////////////////////////////////
        // pz = pT * tanLambda = a * |Bz| / omega * tanLambda //
        ////////////////////////////////////////////////////////

        J[5][0] =  0.0;  
        // dpz / dd0

        J[5][1] =  0.0;  
        // dpz / dphi0

        J[5][2] =  - pz / std::abs(omega);
        // dpz / domega

        J[5][3] =  0.0;  
        // dpz / dz0

        J[5][4] = pt;
        // dpz / dtanLambda

        std::array<double, 36> C_cart;
        C_cart.fill(0.0);
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                double sum = 0.0;
                for (int k = 0; k < 5; ++k) {
                    for (int l = 0; l < 5; ++l) {
                    
                        sum += J[i][k] * C_helix[k*5 + l] * J[j][l];
                    }
                }
                C_cart[i*6 + j] = sum;
            }
        }

        return C_cart;

    }

    /**
    * @brief Computes the point of closest approach (PCA) and related parameters 
    *        of a charged particle's trajectory in a uniform magnetic field.
    *
    * This function approximates the helical motion of a charged particle in a 
    * uniform magnetic field oriented along the z-axis (Bz) and calculates:
    *  - The point on the helix closest to a reference point (`refPoint`) in 3D space.
    *  - The azimuthal angle (`Phi0`) of the trajectory at the PCA.
    *
    * The procedure is as follows:
    * 1. Compute the transverse momentum `pt` and check it is non-zero.
    * 2. Calculate the radius `R` of the helix in the transverse plane using 
    *    R = pt / (0.3 * |q| * Bz) * 100 (pt in GeV/c, Bz in Tesla, R in cm).
    * 3. Determine the center of the circular projection of the helix in the XY plane.
    * 4. Find the unit vector pointing from the circle center to the reference point.
    * 5. Compute the closest point on the circle to the reference point in XY.
    * 6. Compute the tangent vector at the closest point in the XY plane.
    * 7. Calculate the azimuthal angle `Phi0` of the tangent.
    * 8. Compute the z-coordinate of the PCA along the helix using the straight-line approximation.
    *
    * @param position  Initial 3D position of the particle.
    * @param momentum  3D momentum vector of the particle.
    * @param charge    Particle charge (in e units).
    * @param refPoint  Reference point to which the PCA is calculated.
    * @param Bz        Magnetic field along the z-axis (Tesla).
    *
    * @return A `PCAInfoHelper` structure containing:
    *         - `PCA`   : The 3D position of the closest approach point.
    *         - `Phi0`  : Azimuthal angle of the trajectory at the PCA (radians).
    *
    * @throws std::runtime_error if the transverse momentum is zero or the reference point
    *         coincides with the circle center in the XY plane.
    *
    * @note Assumes a uniform magnetic field along z and a simple helical trajectory.
    */
    GenfitTrack::PCAInfoHelper GenfitTrack::PCAInfo(TVector3 position, TVector3 momentum, int charge, TVector3 refPoint, double Bz) {
            
        double px = momentum.X();
        double py = momentum.Y();
        double pt = momentum.Perp();

        if (pt == 0.0) {
            throw std::runtime_error("Transverse momentum is zero");
        }
  
        double R = pt / (0.3 *std::abs(charge) * Bz) * 100; // in cm (assuming pt in GeV/c, Bz in Tesla)

        double tx = px / pt;
        double ty = py / pt;

        double nx = charge * (ty);
        double ny = charge * (-tx);

        double xc = position[0] + R * nx;   
        double yc = position[1] + R * ny;
        TVector3 center(xc, yc, 0.0);

        TVector3 v = refPoint - center;

        double vxy = std::sqrt(v.X()*v.X() + v.Y()*v.Y());

        if (vxy == 0.0) {
            throw std::runtime_error("Reference point coincides with circle center");
        }
            
        TVector3 u(v.X()/vxy, v.Y()/vxy, 0.0);
        TVector3 closestPoint = center + R * u;
        TVector3 r = closestPoint - center;

        int sign = (charge > 0) ? 1 : -1;
        TVector3 tangent(
            -sign * r.Y(),
            sign * r.X(),
            0.0
        );

        tangent = tangent.Unit();

        double angleToX = std::atan2(tangent.Y(), tangent.X());

        double pR   = momentum.Perp();   
        double pZ   = momentum.Z();     
        double R0 = position.Perp();    
        double Z0 = position.Z();

        double tPCA = -(R0*pR + Z0*pZ)/(pR*pR + pZ*pZ);
        double ZPCA = Z0 + pZ * tPCA;

        PCAInfoHelper helper;

        closestPoint.SetZ(ZPCA);
        helper.PCA = closestPoint;

        helper.Phi0 = angleToX;
            
        return helper;

    }

}