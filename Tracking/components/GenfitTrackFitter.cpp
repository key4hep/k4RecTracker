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

//=== Standard Library ===
#include <algorithm>
#include <cmath>
#include <cstdlib>  
#include <fstream>  
#include <filesystem>  
#include <iostream>
#include <iterator> 
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

//=== ROOT / C++ ABI ===
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TFile.h"
#include <cxxabi.h>
#include <Eigen/Dense>

//=== GenFit ===
#include <ConstField.h>
#include <DAF.h>
#include <EventDisplay.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitterRefTrack.h>
#include <MaterialEffects.h>
#include <PlanarMeasurement.h>
#include <RKTrackRep.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>
#include <TGeoMaterialInterface.h>

//=== Gaudi / k4FWCore ===
#include "Gaudi/Property.h"
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"

//=== DD4hep / DDRec / DDSegmentation ===
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DetectorSelector.h"
#include "DD4hep/Readout.h"
#include <DDRec/DetectorData.h>
#include <DDRec/Vector3D.h>
#include <DDSegmentation/BitFieldCoder.h>
#include "DD4hep/Fields.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/MaterialManager.h"

//=== podio / edm4hep ===
#include "podio/Frame.h"
#include "podio/ROOTReader.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackMCParticleLinkCollection.h"

//=== k4Interface ===
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

//=== GenfitInterface ===
#include "GenfitWireMeasurement.h"
#include "GenfitPlanarMeasurement.h"
#include "GenfitField.h"
#include "GenfitTrack.h"
#include "GenfitMaterialInterface.h"

//=== Others ===
#include <marlinutil/HelixClass_double.h>
#include <Objects/Helix.h>

// Define collection types
#include "podio/UserDataCollection.h"

#include "utils.h"
#include "FastCircleSeed.h"



/** @struct GenfitTrackFitter
 *
 *  Gaudi MultiTransformer that refines the parameters of the reconstructed tracks using the GENFIT library.  
 *  The fitting process minimizes the track χ² while accounting for material effects, magnetic field, and detector geometry.
 *
 *  If the fit fails for a given hypothesis, a fallback track is still generated with χ² = -1 and ndf = -1 to preserve collection integrity.
 *  Each fitted track also includes the extrapolated state at the calorimeter (barrel or endcap).
 *
 *  input:
 *    - initial track collection : edm4hep::TrackCollection
 *
 *  output:
 *    - fitted track collection : edm4hep::TrackCollection
 *
 *  @author Andrea De Vita
 *  @date   2026-02
 *
*/

/*
* Debug levels:
* - 0 no debug
* - 1 show trackStates
* - 2 show measurements and trackStates
* - 3 show measurements, trackStates and debugLevel=1 fitter
* - 4 show measurements, trackStates, debugLevel=1 fitter and materialEffects = 1
*
*/

struct GenfitTrackFitter final : 
        k4FWCore::MultiTransformer< std::tuple< edm4hep::TrackCollection>(const edm4hep::TrackCollection&)>                                                                         
{
    GenfitTrackFitter(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                
                KeyValues("InputTracks", {"InputTracks"})
                
            },
            {   

                KeyValues("OutputFittedTracks", {"Fitted_tracks"})

            }) {}
    
            
    StatusCode initialize() {     
        
        if (!gGeoManager) {
            std::cerr << "Error: TGeoManager is not initialized!" << std::endl;
            return StatusCode::FAILURE;
        }
        
        // Initialize the Genfit: FieldManager and MaterialEffects
        m_detector = m_geoSvc->getDetector();
        m_field = m_detector->field();
        m_genfitField=new GenfitInterface::GenfitField(m_field);

        m_fieldManager = genfit::FieldManager::getInstance();
        m_fieldManager->init(m_genfitField);

        m_geoMaterial=GenfitInterface::GenfitMaterialInterface::getInstance(m_detector);      
        // genfit::MaterialEffects::getInstance()->setEnergyLossBrems(true);
        // genfit::MaterialEffects::getInstance()->setNoiseBrems(true);
        // genfit::MaterialEffects::getInstance()->setMscModel("Highland");

        genfit::MaterialEffects::getInstance()->setEnergyLossBrems(false);
        genfit::MaterialEffects::getInstance()->setNoiseBrems(false);   

        int debug_lvl_material = 0;
        if (m_debug_lvl >= 4)
        {
            debug_lvl_material = 1;
        }
        genfit::MaterialEffects::getInstance()->setDebugLvl(debug_lvl_material);         

        // Retrieve the SurfaceManager, ddd4hep::rec::DCH_info and dd4hep::DDSegmentation::BitFieldCoder
        // These object are necessary to extract the drift chamber hits information, such as positions of the wire extremities
        m_surfMan = m_geoSvc->getDetector()->extension<dd4hep::rec::SurfaceManager>();
        // If the detector doesn't have a drift chamber, this part will be skipped
        try {
            std::string DCH_name("DCH_v2");
            dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
            m_dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();
            dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector("DCH_v2");
            dd4hep::Readout dch_readout = dch_sd.readout();
            m_dc_decoder = dch_readout.idSpec().decoder();
        } catch (const std::out_of_range& e) {}

        // Retrive calorimeter information
        // These parameters are necessary to propagate the track to the calorimeter surface
        dd4hep::Detector& detector = dd4hep::Detector::getInstance();
        const std::vector< dd4hep::DetElement>& isDuaReadout = dd4hep::DetectorSelector(detector).detectors(  dd4hep::DetType::CALORIMETER | dd4hep::DetType::BARREL | dd4hep::DetType::ENDCAP, dd4hep::DetType::AUXILIARY );
        if( isDuaReadout.size() > 0 )
        {

            const dd4hep::rec::LayeredCalorimeterData * DualReadoutExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));

            // mm (dd4hep::mm = 0.1)
            m_eCalBarrelInnerR  = DualReadoutExtension->extent[0] / dd4hep::mm;         // barrel rmin
            m_eCalBarrelMaxZ    = DualReadoutExtension->extent[2] / dd4hep::mm;         // barrel zmax == endcap zmin

            m_eCalEndCapInnerR  = DualReadoutExtension->extent[4] / dd4hep::mm;         // endcap rmin
            m_eCalEndCapOuterR  = DualReadoutExtension->extent[5] / dd4hep::mm;         // endcap rmax
            m_eCalEndCapInnerZ  = DualReadoutExtension->extent[2] / dd4hep::mm;         // endcap zmin
            m_eCalEndCapOuterZ  = DualReadoutExtension->extent[3] / dd4hep::mm;         // endcap zmax
        }
        else
        {
            const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
                                                                                              ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));

            // mm (dd4hep::mm = 0.1)                                                                                  
            m_eCalBarrelInnerR = eCalBarrelExtension->extent[0] / dd4hep::mm;
            m_eCalBarrelMaxZ = eCalBarrelExtension->extent[3] / dd4hep::mm;
                    
            const dd4hep::rec::LayeredCalorimeterData * eCalEndCapExtension = getExtension( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
                                                                                              ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ));


            // mm (dd4hep::mm = 0.1)                                                                                  
            m_eCalEndCapInnerR = eCalEndCapExtension->extent[0] / dd4hep::mm;
            m_eCalEndCapOuterR = eCalEndCapExtension->extent[1] / dd4hep::mm;
            m_eCalEndCapInnerZ = eCalEndCapExtension->extent[2] / dd4hep::mm;
            m_eCalEndCapOuterZ = eCalEndCapExtension->extent[3] / dd4hep::mm;
        }

        // N.B. we are assuming that the magnetic field is uniform and along the z direction
        // Z-component of the magnetic field at the center of the detector
        // This component is used to propagate the track to the calorimeter surface, but it is also necessary to
        // compute the pt component given omega
        dd4hep::Position center(0, 0, 0);
        dd4hep::Direction bfield = m_field.magneticField(center);
        m_Bz = bfield.z() / dd4hep::tesla;

        return StatusCode::SUCCESS;

    }
    
    std::tuple<edm4hep::TrackCollection> operator()( const edm4hep::TrackCollection& tracks_input) const override                                                 
    {
        
        // This collection stores the output of the fit
        edm4hep::TrackCollection FittedTracks;

        info() << "Event number: " << event_counter++ << endmsg;

        // Loop over the tracks created by the pattern recognition step
        for (const auto& track : tracks_input)
        {

            num_tracks  +=1;

            // Skip background tracks if the option is enabled
            // Consider background tracks those with type = 0
            if (m_skip_background && track.getType() == 0) 
            {
                num_skip += 1;
                warning() << "Track " << num_tracks - 1 << ": background track (type = 0), skipping fit.\n" << endmsg;
                continue;        // skip background        
            } 
            
            if (track.getTrackerHits().size() <= 3) 
            {   
                num_skip += 1;
                warning() << "Track " << num_tracks - 1 << ": less than 3 hits, skipping fit.\n" << endmsg;
                continue;        // skip tracks with less then 3 hits (seed initialization needs 3 hits)
            }

            num_processed_tracks += 1;

            int isSuccess = 0;
            if (m_singleEvaluation)
            {

                
                GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(track, m_dch_info, m_dc_decoder);
                track_interface.InitializeTrack(m_Bz, false, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt);  
                
                auto track_init = track_interface.GetInitialization();
                if (m_debug_lvl > 0)
                {
                    debug() << "Track " << num_tracks - 1 << " with " << track.getTrackerHits().size()<< " hits: initial seed for track fit:" << endmsg;
                    debug() << "  Initial position [mm]: (" << track_init.Position.X() / dd4hep::mm << ", " << track_init.Position.Y() / dd4hep::mm << ", " << track_init.Position.Z() / dd4hep::mm << ")" << endmsg;
                    debug() << "  Initial momentum [GeV/c]: (" << track_init.Momentum.X() << ", " << track_init.Momentum.Y() << ", " << track_init.Momentum.Z() << ")" << endmsg;
                    debug() << "  Initial Covariance Matrix:" << endmsg;
                    track_init.CovMatrix.Print();
                    debug() << "  Charge hypothesis: " << track_init.Charge << endmsg;
                    debug() << "  Max hit for loopers: " << track_init.NumHits << endmsg;
                    debug() << "\n" << endmsg;
                }

                int debug_track = 0;
                if (m_debug_lvl > 1) debug_track = 1;
                int pdgCode = m_particleHypotesis[0];    
                track_interface.CreateGenFitTrack(pdgCode, debug_track);   
    

                bool isFit = track_interface.Fit(m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, m_debug_lvl);  

                if (isFit)
                {

                    isSuccess = 1;

                    auto edm4hep_track = track_interface.GetTrack_edm4hep();

                    // Propagate to calorimeter and store the state at calorimeter
                    FillTrackWithCalorimeterExtrapolation(  edm4hep_track, m_Bz, track_interface.GetCharge(), 
                                                            m_eCalBarrelInnerR, m_eCalBarrelMaxZ, m_eCalEndCapInnerR, m_eCalEndCapOuterR, m_eCalEndCapInnerZ);
                    
                    
                    if (m_debug_lvl > 0)
                    {   
                        auto trackStates = edm4hep_track.getTrackStates();
                        edm4hep::TrackState trackStateCalo;
                        for (const auto& ts : trackStates)
                        {
                            if (ts.location == edm4hep::TrackState::AtCalorimeter)
                            {
                                trackStateCalo = ts;
                                break;
                            }
                        }

                        std::cout << "GenfitTrackFitter    DEBUG : TrackState at Calo: " << std::endl;
                        std::cout << "  D0: " << trackStateCalo.D0 << " mm" << std::endl;
                        std::cout << "  Z0: " << trackStateCalo.Z0 << " mm" << std::endl;
                        std::cout << "  phi: " << trackStateCalo.phi << " rad" << std::endl;
                        std::cout << "  omega: " << trackStateCalo.omega << " 1/mm" << std::endl;
                        std::cout << "  tanLambda: " << trackStateCalo.tanLambda << std::endl;
                        std::cout << "  location: " << trackStateCalo.location << std::endl;
                        std::cout << "  reference point: (" << trackStateCalo.referencePoint.x << ", " << trackStateCalo.referencePoint.y << ", " << trackStateCalo.referencePoint.z << ") mm" << std::endl;
                        
                    }

                    // Add the fitted track to the output collection
                    FittedTracks.push_back(edm4hep_track);

                }
                

            }
            else
            {

                int winning_hypothesis = -1; double winning_chi2_ndf = std::numeric_limits<double>::max();
                for (int pdgCode : m_particleHypotesis)
                {

                    GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(track, m_dch_info, m_dc_decoder);
                    track_interface.InitializeTrack(m_Bz, false, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt);   

                    track_interface.CreateGenFitTrack(pdgCode, 0);   

                    bool isFit = track_interface.Fit(m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, 0);  

                    if (isFit)
                    {

                        auto edm4hep_track = track_interface.GetTrack_edm4hep();
                        float genfit_chi2_val = edm4hep_track.getChi2();
                        int genfit_ndf_val = edm4hep_track.getNdf();

                        if (genfit_chi2_val <= 0 || genfit_ndf_val <= 0)  continue; // skip invalid fits
                        
                        double chi2_ndf = genfit_chi2_val / genfit_ndf_val;
                        if (chi2_ndf < winning_chi2_ndf)
                        {
                            winning_chi2_ndf = chi2_ndf;
                            winning_hypothesis = pdgCode;
                        }
                    
                    }

                }

                if (winning_hypothesis == -1)
                {
                    debug() << "Track " << num_tracks - 1 << ": fit failed for all hypotheses, trying with less hits." << endmsg;
                }
                else
                {

                    isSuccess = 1;
                    debug() << "Track " << num_tracks - 1 << ": winning hypothesis is " << winning_hypothesis << " with chi2 / ndf = " << winning_chi2_ndf << endmsg;
                    int pdgCode = winning_hypothesis;

                    GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(track, m_dch_info, m_dc_decoder);
           
                    track_interface.InitializeTrack(m_Bz, false, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt);   

                    int debug_track = 0;
                    if (m_debug_lvl > 1) debug_track = 1;
                    track_interface.CreateGenFitTrack(pdgCode, debug_track);   
        
                    bool isFit = track_interface.Fit(m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, m_debug_lvl);

                    // Extra check
                    if (!isFit)
                    {
                        std::cerr << "Error: track fitting failed unexpectedly." << std::endl;
                        std::exit(EXIT_FAILURE);
                    }

                    auto edm4hep_track = track_interface.GetTrack_edm4hep();

                    // Propagate to calorimeter and store the state at calorimeter
                    FillTrackWithCalorimeterExtrapolation(  edm4hep_track, m_Bz, track_interface.GetCharge(), 
                                                            m_eCalBarrelInnerR, m_eCalBarrelMaxZ, m_eCalEndCapInnerR, m_eCalEndCapOuterR, m_eCalEndCapInnerZ);
                    
                    if (m_debug_lvl > 0)
                    {   
                        auto trackStates = edm4hep_track.getTrackStates();
                        edm4hep::TrackState trackStateCalo;
                        for (const auto& ts : trackStates)
                        {
                            if (ts.location == edm4hep::TrackState::AtCalorimeter)
                            {
                                trackStateCalo = ts;
                                break;
                            }
                        }

                        std::cout << "GenfitTrackFitter    DEBUG : TrackState at Calo: " << std::endl;
                        std::cout << "  D0: " << trackStateCalo.D0 << " mm" << std::endl;
                        std::cout << "  Z0: " << trackStateCalo.Z0 << " mm" << std::endl;
                        std::cout << "  phi: " << trackStateCalo.phi << " rad" << std::endl;
                        std::cout << "  omega: " << trackStateCalo.omega << " a.u." << std::endl;
                        std::cout << "  tanLambda: " << trackStateCalo.tanLambda << std::endl;
                        std::cout << "  location: " << trackStateCalo.location << std::endl;
                        std::cout << "  reference point: (" << trackStateCalo.referencePoint.x << ", " << trackStateCalo.referencePoint.y << ", " << trackStateCalo.referencePoint.z << ") mm" << std::endl;
                            
                    }
                    
                    // Add the fitted track to the output collection
                    FittedTracks.push_back(edm4hep_track);

                }

            }

            if (!isSuccess)
            {
                if (m_singleEvaluation)
                {
                    GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(track, m_dch_info, m_dc_decoder);
                    track_interface.InitializeTrack(m_Bz, true, 1, std::nullopt, std::nullopt, 0.05, 3);  
                
                    auto track_init = track_interface.GetInitialization();
                    if (m_debug_lvl > 0)
                    {
                        debug() << "Track " << num_tracks - 1 << " with " << track.getTrackerHits().size() << " hits: initial seed for track fit:" << endmsg;
                        debug() << "  Initial position [mm]: (" << track_init.Position.X() / dd4hep::mm << ", " << track_init.Position.Y() / dd4hep::mm << ", " << track_init.Position.Z() / dd4hep::mm << ")" << endmsg;
                        debug() << "  Initial momentum [GeV/c]: (" << track_init.Momentum.X() << ", " << track_init.Momentum.Y() << ", " << track_init.Momentum.Z() << ")" << endmsg;
                        debug() << "  Initial Covariance Matrix:" << endmsg;
                        track_init.CovMatrix.Print();
                        debug() << "  Charge hypothesis: " << track_init.Charge << endmsg;
                        debug() << "  Max hit for loopers: " << track_init.NumHits << endmsg;
                        debug() << "\n" << endmsg;
                    }

                    int debug_track = 0;
                    if (m_debug_lvl > 1) debug_track = 1;
                    int pdgCode = m_particleHypotesis[0];    
                    track_interface.CreateGenFitTrack(pdgCode, debug_track);   
        

                    bool isFit = track_interface.Fit(m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, m_debug_lvl);  

                    if (isFit)
                    {

                        isSuccess = 1;

                        auto edm4hep_track = track_interface.GetTrack_edm4hep();

                        // Propagate to calorimeter and store the state at calorimeter
                        FillTrackWithCalorimeterExtrapolation(  edm4hep_track, m_Bz, track_interface.GetCharge(), 
                                                                m_eCalBarrelInnerR, m_eCalBarrelMaxZ, m_eCalEndCapInnerR, m_eCalEndCapOuterR, m_eCalEndCapInnerZ);
                        
                        
                        if (m_debug_lvl > 0)
                        {   
                            auto trackStates = edm4hep_track.getTrackStates();
                            edm4hep::TrackState trackStateCalo;
                            for (const auto& ts : trackStates)
                            {
                                if (ts.location == edm4hep::TrackState::AtCalorimeter)
                                {
                                    trackStateCalo = ts;
                                    break;
                                }
                            }

                            std::cout << "GenfitTrackFitter    DEBUG : TrackState at Calo: " << std::endl;
                            std::cout << "  D0: " << trackStateCalo.D0 << " mm" << std::endl;
                            std::cout << "  Z0: " << trackStateCalo.Z0 << " mm" << std::endl;
                            std::cout << "  phi: " << trackStateCalo.phi << " rad" << std::endl;
                            std::cout << "  omega: " << trackStateCalo.omega << " 1/mm" << std::endl;
                            std::cout << "  tanLambda: " << trackStateCalo.tanLambda << std::endl;
                            std::cout << "  location: " << trackStateCalo.location << std::endl;
                            std::cout << "  reference point: (" << trackStateCalo.referencePoint.x << ", " << trackStateCalo.referencePoint.y << ", " << trackStateCalo.referencePoint.z << ") mm" << std::endl;
                            
                        }

                        // Add the fitted track to the output collection
                        FittedTracks.push_back(edm4hep_track);


                    }
                    else
                    {

                        number_failures += 1;
            
                        debug() << "Track " << num_tracks - 1 << ": fit failed for single evaluation hypothesis, skipping track." << endmsg;
                        auto failedTrack = FittedTracks.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                        continue;
                    }
        
                    
                }
                else
                {

                    int winning_hypothesis = -1; double winning_chi2_ndf = std::numeric_limits<double>::max();
                    for (int pdgCode : m_particleHypotesis)
                    {

                        GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(track, m_dch_info, m_dc_decoder);
                        track_interface.InitializeTrack(m_Bz, true, 1, std::nullopt, std::nullopt, 0.05, 3);   

                        track_interface.CreateGenFitTrack(pdgCode, 0);   

                        bool isFit = track_interface.Fit(m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, 0);  

                        if (isFit)
                        {

                            auto edm4hep_track = track_interface.GetTrack_edm4hep();
                            float genfit_chi2_val = edm4hep_track.getChi2();
                            int genfit_ndf_val = edm4hep_track.getNdf();

                            if (genfit_chi2_val <= 0 || genfit_ndf_val <= 0)  continue; // skip invalid fits
                            
                            double chi2_ndf = genfit_chi2_val / genfit_ndf_val;
                            if (chi2_ndf < winning_chi2_ndf)
                            {
                                winning_chi2_ndf = chi2_ndf;
                                winning_hypothesis = pdgCode;
                            }
                        
                        }

                    }

                    if (winning_hypothesis == -1)
                    {

                        debug() << "Track " << num_tracks - 1 << ": fit failed for all hypotheses." << endmsg;
                        number_failures += 1;
                        auto failedTrack = FittedTracks.create();
                        failedTrack.setChi2(-1);
                        failedTrack.setNdf(-1);
                        continue;

                    }
                    else
                    {

                        debug() << "Track " << num_tracks - 1 << ": winning hypothesis is " << winning_hypothesis << " with chi2 / ndf = " << winning_chi2_ndf << endmsg;
                        int pdgCode = winning_hypothesis;

                        GenfitInterface::GenfitTrack track_interface = GenfitInterface::GenfitTrack(track, m_dch_info, m_dc_decoder);
            
                        track_interface.InitializeTrack(m_Bz, false, 1, std::nullopt, std::nullopt, std::nullopt, std::nullopt);   

                        int debug_track = 0;
                        if (m_debug_lvl > 1) debug_track = 1;
                        track_interface.CreateGenFitTrack(pdgCode, debug_track);   
            

                        bool isFit = track_interface.Fit(m_Beta_init, m_Beta_final, m_Beta_steps, m_Bz, m_debug_lvl);

                        // Extra check
                        if (!isFit)
                        {
                            std::cerr << "Error: track fitting failed unexpectedly." << std::endl;
                            std::exit(EXIT_FAILURE);
                        }

                        auto edm4hep_track = track_interface.GetTrack_edm4hep();

                        // Propagate to calorimeter and store the state at calorimeter
                        FillTrackWithCalorimeterExtrapolation(  edm4hep_track, m_Bz, track_interface.GetCharge(), 
                                                                m_eCalBarrelInnerR, m_eCalBarrelMaxZ, m_eCalEndCapInnerR, m_eCalEndCapOuterR, m_eCalEndCapInnerZ);
                        
                        
                        if (m_debug_lvl > 0)
                        {   
                            auto trackStates = edm4hep_track.getTrackStates();
                            edm4hep::TrackState trackStateCalo;
                            for (const auto& ts : trackStates)
                            {
                                if (ts.location == edm4hep::TrackState::AtCalorimeter)
                                {
                                    trackStateCalo = ts;
                                    break;
                                }
                            }

                            std::cout << "GenfitTrackFitter    DEBUG : TrackState at Calo: " << std::endl;
                            std::cout << "  D0: " << trackStateCalo.D0 << " mm" << std::endl;
                            std::cout << "  Z0: " << trackStateCalo.Z0 << " mm" << std::endl;
                            std::cout << "  phi: " << trackStateCalo.phi << " rad" << std::endl;
                            std::cout << "  omega: " << trackStateCalo.omega << " a.u." << std::endl;
                            std::cout << "  tanLambda: " << trackStateCalo.tanLambda << std::endl;
                            std::cout << "  location: " << trackStateCalo.location << std::endl;
                            std::cout << "  reference point: (" << trackStateCalo.referencePoint.x << ", " << trackStateCalo.referencePoint.y << ", " << trackStateCalo.referencePoint.z << ") mm" << std::endl;
                                
                        }
                        
                        // Add the fitted track to the output collection
                        FittedTracks.push_back(edm4hep_track);

                    }
                }
            }

        }

        return std::make_tuple( std::move(FittedTracks));
        
    } 
    
    StatusCode finalize() {     
        
        info() << "Run report:" << endmsg;
        info() << "Number of tracks: " << num_tracks << endmsg;
        info() << "Number of successes: " << (num_processed_tracks - number_failures) <<  "/" << num_processed_tracks << endmsg;
        info() << "Number of skipped tracks: " << num_skip << "/" << num_tracks << endmsg;
        info() << "----------------\n" << endmsg;

        return StatusCode::SUCCESS;

    }

    public:
        
        mutable int event_counter = 0;

        // Num_tracks = num_processed_tracks + num_skip
        mutable int num_tracks = 0;             // Total number of tracks
        mutable int num_skip = 0;               // Number of skipped tracks
        mutable int num_processed_tracks = 0;   // Number of tracks that have been processed (not skipped)
        mutable int number_failures = 0;        // Number of failed fits


    private:

        

        ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
        dd4hep::Detector* m_detector{nullptr};  // Detector instance
        dd4hep::OverlayedField m_field;         // Magnetic field

        GenfitInterface::GenfitField* m_genfitField;
        genfit::FieldManager* m_fieldManager;

        GenfitInterface::GenfitMaterialInterface* m_geoMaterial;
        genfit::MaterialEffects* m_materialEffects;

        dd4hep::rec::SurfaceManager* m_surfMan;
        dd4hep::rec::DCH_info* m_dch_info;
        dd4hep::DDSegmentation::BitFieldCoder* m_dc_decoder;

        double m_Bz;

        double m_eCalBarrelInnerR;
        double m_eCalBarrelMaxZ;
        double m_eCalEndCapInnerR;
        double m_eCalEndCapOuterR;
        double m_eCalEndCapInnerZ;
        double m_eCalEndCapOuterZ;

        std::vector<int> m_particleHypotesis = {211};   // {11,13,211,321,2212} -> e, mu, pi, K, p

        Gaudi::Property<double> m_Beta_init{this, "Beta_init", 100, "Beta Initial value"};
        Gaudi::Property<double> m_Beta_final{this, "Beta_final", 0.05, "Beta Final value"};
        Gaudi::Property<int> m_Beta_steps{this, "Beta_steps", 15, "Beta number of Steps"};


        Gaudi::Property<bool> m_skip_background{this, "skip_background", true, "Skip background track"};
        Gaudi::Property<bool> m_singleEvaluation{this, "single_evaluation", false, "Single evaluation mode"};
        


        Gaudi::Property<int> m_debug_lvl{this, "debug_lvl", 0, "Debug level"}; 

};

DECLARE_COMPONENT(GenfitTrackFitter)
