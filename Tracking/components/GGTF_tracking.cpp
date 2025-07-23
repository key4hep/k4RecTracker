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

// Standard Library
#include <algorithm>
#include <cmath>
#include <cstdlib>         // For getenv
#include <filesystem>      // For std::filesystem::path
#include <fstream>         // For std::ifstream
#include <iostream>
#include <iterator>        // For std::istreambuf_iterator
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

// ONNX & Torch
#include <ATen/ATen.h>
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"
#include "onnxruntime_run_options_config_keys.h"

// ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TGeoMatrix.h"
#include "TRandom3.h"
#include "TVector3.h"

// === Gaudi Framework ===
#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// === k4FWCore / k4Interface ===
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// === EDM4HEP & PODIO ===
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "podio/UserDataCollection.h"

// === EDM4HEP Extensions ===
#include "extension/DriftChamberDigiCollection.h"
#include "extension/DriftChamberDigiLocalCollection.h"
#include "extension/MCRecoDriftChamberDigiAssociationCollection.h"
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"
#include "extension/TrackCollection.h"
#include "extension/TrackerHit.h"

// === DD4hep ===
#include "DD4hep/Detector.h"       
#include "DDRec/DCH_info.h"
#include "DDRec/Vector3D.h"
#include "DDSegmentation/BitFieldCoder.h"

// === Project-specific ===
#include "utils.hpp"

/** @struct GGTF_tracking
 *
 *  Gaudi MultiTransformer that generates a Track collection by analyzing the digitalized hits through the GGTF_tracking. 
 *  The first step takes the raw hits and it returns a collection of 4-dimensional points inside an embedding space.
 *  Eeach 4-dim point has 3 geometric coordinates and 1 charge, the meaning of which can be described intuitively by a potential, 
 *  which attracts hits belonging to the same cluster and drives away those that do not.
 *  This collection of 4-dim points is analysed by a clustering step, which groups together hits belonging to the same track.
 *
 *  input: 
 *    - digitalized hits from DC (global coordinates) : extension::SenseWireHitCollection
 *    - digitalized hits from vertex (global coordinates) : edm4hep::TrackerHitPlaneCollection
 *    - digitalized hits from silicon wrapper (global coordinates) : edm4hep::TrackerHitPlaneCollection
 *
 *  output:
 *    - Track collection : extension::TrackCollection
 *
 *
 *
 *  @author Andrea De Vita, Maria Dolores Garcia, Brieuc Francois
 *  @date   2025-06
 *
 */

struct GGTF_tracking final : 
        k4FWCore::MultiTransformer< std::tuple<extension::TrackCollection>(   
                                                                    const std::vector<const edm4hep::TrackerHitPlaneCollection*>&, 
                                                                    const std::vector<const extension::SenseWireHitCollection*>&)>                                                                                            
{
    GGTF_tracking(const std::string& name, ISvcLocator* svcLoc) : 
        MultiTransformer ( name, svcLoc,
            {   
                 
                KeyValues("inputPlanarHits", {"inputPlanarHits"}),
                KeyValues("inputWireHits", {"inputWireHits"})

            },
            {   
               
               KeyValues("outputTracks", {"outputTracks"})      
            
            }) {m_geoSvc = serviceLocator()->service(m_geoSvcName);}
    
    StatusCode initialize() {

        ///////////////////////////////
        ///// ONNX Initialization /////
        ///////////////////////////////

        // Initialize the ONNX memory info object for CPU memory allocation.
        // This specifies that the memory will be allocated using the Arena Allocator on the CPU.
        fInfo = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

        // Create and initialize the ONNX environment with a logging level set to WARNING.
        // This environment handles logging and runtime configuration for the ONNX session.
        auto envLocal = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
        fEnv          = std::move(envLocal);

        // Set the session options to configure the ONNX inference session.
        // Set the number of threads used for intra-op parallelism to 1 (single-threaded execution).
        fSessionOptions.SetIntraOpNumThreads(1);

        // Disable all graph optimizations to keep the model execution as close to the original as possible.
        fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);
        // fSessionOptions.DisableMemPattern();

        // Create an ONNX inference session using the configured environment and session options.
        // The session is used to load the model specified by the `modelPath`.
        auto sessionLocal = std::make_unique<Ort::Session>(*fEnv, modelPath.value().c_str(), fSessionOptions);
        fSession          = std::move(sessionLocal);

        // Create an ONNX allocator with default options to manage memory allocations during runtime.
        Ort::AllocatorWithDefaultOptions allocator;

        // Get the name of the first input node (index i) from the ONNX model and store it.
        // This retrieves the name from the model and releases the memory after storing it in fInames.
        // Get the name of the first output node (index i) from the ONNX model and store it.
        // This retrieves the name from the model and releases the memory after storing it in fOnames.
        std::size_t i = 0;
        const auto input_name = fSession->GetInputNameAllocated(i, allocator).release();
        const auto output_names = fSession->GetOutputNameAllocated(i, allocator).release();

        // Store the retrieved input and output names in the respective vectors.
        fInames.push_back(input_name);
        fOnames.push_back(output_names);

        //////////////////////////////
        ///// Drift Chamber Info /////
        //////////////////////////////

        // Retrieve the drift chamber information and decoder
        try {
            std::string DCH_name("DCH_v2");
            dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);
            dch_info = DCH_DE.extension<dd4hep::rec::DCH_info>();
            dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector("DCH_v2");
            dd4hep::Readout dch_readout = dch_sd.readout();
            dc_decoder = dch_readout.idSpec().decoder();
        } catch (const std::out_of_range& e) {}



        return StatusCode::SUCCESS;

   }

    
    std::tuple<extension::TrackCollection> operator()(      const std::vector<const edm4hep::TrackerHitPlaneCollection*>& inputPlanarHits, 
                                                            const std::vector<const extension::SenseWireHitCollection*>& inputWireHits) const override 
    {

        ////////////////////////////////////////
        ////////// DATA PREPROCESSING //////////
        ////////////////////////////////////////

        info() << "Processing event number: " << index_counter++ << endmsg;

        // Vector to store the global input values for all hits.
        // This will contain position and other hit-specific data to be used as input for the model.
        std::vector<float> ListGlobalInputs; int globalHitIndex = 0;

        /// Processing hits from the Silicon Detectors
        std::vector<int64_t> ListHitType_Planar; int planarHit_index = 0;
        std::vector<int64_t> ListPlanarHitIndices; int planarHitCollection_index = 0; 
        for (const auto& inputHits_collection : inputPlanarHits)
        {   
            
            int planarHit_subcollection_index = 0;
            for (const auto& input_hit : *inputHits_collection) {
                
                

                // Add the 3D position of the hit to the global input list.
                ListGlobalInputs.push_back(input_hit.getPosition().x);
                ListGlobalInputs.push_back(input_hit.getPosition().y);
                ListGlobalInputs.push_back(input_hit.getPosition().z);
                
                // Add placeholder values for additional input dimensions.
                ListGlobalInputs.push_back(1.0); 
                ListGlobalInputs.push_back(0.0);
                ListGlobalInputs.push_back(0.0);
                ListGlobalInputs.push_back(0.0); 
                
                // Store the current index in ListHitType_Planar and increment the global iterator.
                ListHitType_Planar.push_back(globalHitIndex);

                ListPlanarHitIndices.push_back(planarHitCollection_index);
                ListPlanarHitIndices.push_back(planarHit_subcollection_index);

                globalHitIndex += 1;  
                planarHit_index += 1;   
                planarHit_subcollection_index += 1;   
                
                
            }

            planarHitCollection_index += 1;
        }
        // Convert ListHitType_Planar to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_Planar_tensor = torch::from_blob(ListHitType_Planar.data(), {planarHit_index}, torch::kFloat32);

        int64_t n_planar_hits = ListPlanarHitIndices.size() / 2;
        torch::Tensor ListPlanarHitIndices_tensor = torch::from_blob(
            ListPlanarHitIndices.data(),
            {n_planar_hits, 2},
            torch::kInt64
        );

        /// Processing hits from gaseous detectors
        std::vector<int64_t> ListHitType_Wire; int wireHit_index = 0;
        std::vector<int64_t> ListWireHitIndices; int wireHitCollection_index = 0; 
        for (const auto& inputHits_collection : inputWireHits)
        {
            int wireHit_subcollection_index = 0;
            for (const auto& input_hit : *inputHits_collection) {
                
                // position along the wire, wire direction and drift distance
                edm4hep::Vector3d wirePos = input_hit.getPosition();   
                TVector3 wire_pos(static_cast<float>(wirePos.x),static_cast<float>(wirePos.y),static_cast<float>(wirePos.z));

                double distanceToWire = input_hit.getDistanceToWire(); 
                double wireAzimuthalAngle = input_hit.getWireAzimuthalAngle();  
                double wireStereoAngle = input_hit.getWireStereoAngle();  
                
                // Direction of the wire: z'
                TVector3 direction(0,0,1);
                direction.RotateX(wireStereoAngle);
                direction.RotateZ(wireAzimuthalAngle);
            
                TVector3 z_prime;
                z_prime = direction.Unit();

                // x' axis, orthogonal to z'
                TVector3 x_prime(1.0, 0.0, -direction.X() / direction.Z());
                x_prime = x_prime.Unit();

                // y' axis = z' × x'
                TVector3 y_prime = z_prime.Cross(x_prime);
                y_prime = y_prime.Unit();

                // Define the local left/right positions
                TVector3 left_local_pos(-distanceToWire, 0.0, 0.0);
                TVector3 right_local_pos(distanceToWire, 0.0, 0.0);

                // Transform to global
                TVector3 left_global_pos = x_prime * left_local_pos.X() + y_prime * left_local_pos.Y() + z_prime * left_local_pos.Z() + wire_pos;
                TVector3 right_global_pos = x_prime * right_local_pos.X() + y_prime * right_local_pos.Y() + z_prime * right_local_pos.Z() + wire_pos;
        
                // Add the 3D position of the left hit to the global input list.
                ListGlobalInputs.push_back(left_global_pos.X());
                ListGlobalInputs.push_back(left_global_pos.Y());
                ListGlobalInputs.push_back(left_global_pos.Z());
                
                // Add the difference between the right and left hit positions to the global input list.
                ListGlobalInputs.push_back(0.0); 
                ListGlobalInputs.push_back(right_global_pos.X()-left_global_pos.X());
                ListGlobalInputs.push_back(right_global_pos.y()-left_global_pos.Y());
                ListGlobalInputs.push_back(right_global_pos.Z()-left_global_pos.Z());
                
                // Store the current index in ListHitType_Wire and increment the global iterator.
                ListHitType_Wire.push_back(globalHitIndex);

                ListWireHitIndices.push_back(wireHitCollection_index);
                ListWireHitIndices.push_back(wireHit_subcollection_index);

                globalHitIndex += 1;  
                wireHit_index += 1;   
                wireHit_subcollection_index += 1;                  
            }
        }
        // Convert ListHitType_CDC to a Torch tensor for use in PyTorch models.
        torch::Tensor ListHitType_Wire_tensor = torch::from_blob(ListHitType_Wire.data(), {wireHit_index}, torch::kFloat32);

        int64_t n_wire_hits = ListWireHitIndices.size() / 2;
        torch::Tensor ListWireHitIndices_tensor = torch::from_blob(
            ListWireHitIndices.data(),
            {n_wire_hits, 2},
            torch::kInt64
        );

        torch::Tensor planar_type = torch::zeros({n_planar_hits, 1}, torch::kInt64);
        torch::Tensor wire_type = torch::ones({n_wire_hits, 1}, torch::kInt64);

        torch::Tensor planar_with_type = torch::cat({planar_type, ListPlanarHitIndices_tensor}, 1);
        torch::Tensor wire_with_type   = torch::cat({wire_type,   ListWireHitIndices_tensor}, 1);

        torch::Tensor ListHitIndicesGlobal = torch::cat({planar_with_type, wire_with_type}, 0);


        //////////////////////////////////
        ////////// TRACK FINDER //////////
        //////////////////////////////////

        // Create a new TrackCollection and TrackerHit3DCollection for storing the output tracks and hits
        extension::TrackCollection output_tracks;
        constexpr int kMaxHits = 20000;
        if (globalHitIndex > 0 && globalHitIndex < kMaxHits)
        {

            ///////////////////////////////
            ////////// GATR STEP //////////
            ///////////////////////////////

            // Calculate the total size of the input tensor, based on the number of hits (it) and the 
            // number of features per hit (7: x, y, z, and four placeholders).
            size_t total_size = globalHitIndex * 7;
            std::vector<int64_t> tensor_shape = {globalHitIndex, 7};

            // Create a vector to store the input tensors that will be fed into the ONNX model.
            std::vector<Ort::Value> input_tensors;
            input_tensors.emplace_back(Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));
            
            // Run the ONNX inference session with the provided input tensor.
            auto output_model_tensors = fSession->Run(Ort::RunOptions{nullptr}, fInames.data(), input_tensors.data(), fInames.size(), fOnames.data(), fOnames.size());
            float* floatarr = output_model_tensors.front().GetTensorMutableData<float>();
            std::vector<float> output_model_vector(floatarr, floatarr + globalHitIndex * 4);

            /////////////////////////////////////
            ////////// CLUSTERING STEP //////////
            /////////////////////////////////////

            auto clustering = get_clustering(output_model_vector, globalHitIndex,tbeta,td);
            torch::Tensor unique_tensor;
            torch::Tensor inverse_indices;
            std::tie(unique_tensor, inverse_indices) = at::_unique(clustering, true, true);
          
            /////////////////////////////////
            ////////// OUTPUT STEP //////////
            /////////////////////////////////

            // Get the total number of unique tracks based on the unique_tensor size
            int64_t number_of_tracks = unique_tensor.numel(); 
            bool has_zero = (unique_tensor == 0).any().item<bool>();
            if (!has_zero)
            {
                auto output_track = output_tracks.create();
                output_track.setType(0);
            }

            // Loop through each unique track ID
            for (int i = 0; i < number_of_tracks; ++i) {

                // Retrieve the current track ID
                auto id_of_track = unique_tensor.index({i});


                // Create a new track in the output collection and set its type to the current track ID
                auto output_track = output_tracks.create();
                output_track.setType(id_of_track.item<int>());

                // Create a mask to select all hits belonging to the current track
                torch::Tensor mask = (clustering == id_of_track);
                
                // Find the indices of the hits that belong to the current track
                torch::Tensor indices = torch::nonzero(mask).flatten();

                auto ListHitIndicesGlobal_view = ListHitIndicesGlobal.accessor<int64_t, 2>();
                auto indices_view  = indices.accessor<int64_t, 1>();

                for (int64_t j = 0; j < indices.size(0); ++j) {
                    int64_t row = indices_view[j];
                    int64_t type = ListHitIndicesGlobal_view[row][0];
                    int64_t idx1 = ListHitIndicesGlobal_view[row][1];
                    int64_t idx2 = ListHitIndicesGlobal_view[row][2];

                    if (type == 0) {
                        // planar hit
                        auto planar_collection = inputPlanarHits[idx1];
                        auto hit = planar_collection->at(idx2);
                        output_track.addToTrackerHits(hit);
                        
                    } else {
                        // wire hit
                        auto wire_collection = inputWireHits[idx1];
                        auto hit = wire_collection->at(idx2);
                        output_track.addToTrackerHits(hit);
                        
                    }
                }

            }

            inverse_indices.reset();
            unique_tensor.reset();
            clustering.reset();
            
            input_tensors.clear();
            output_model_tensors.clear();
            
            
        }

        ListHitType_Planar_tensor.reset();
        ListPlanarHitIndices_tensor.reset();
        ListHitType_Wire_tensor.reset();
        ListWireHitIndices_tensor.reset();
        

        std::vector<int64_t>().swap(ListHitType_Planar);
        std::vector<int64_t>().swap(ListPlanarHitIndices);
        std::vector<int64_t>().swap(ListHitType_Wire);
        std::vector<int64_t>().swap(ListWireHitIndices);

        // Return the output collections as a tuple
        return std::make_tuple(std::move(output_tracks));


    } 

    StatusCode finalize() {     
        
        info() << "Run report:" << endmsg;
        info() << "Number of analysed events: " << index_counter << endmsg;
        info() << "----------------\n" << endmsg;



        return StatusCode::SUCCESS;

    }

    public:
        mutable int index_counter = 0;

    private:

        /// Pointer to the ONNX environment.
        /// This object manages the global state of the ONNX runtime, such as logging and threading.
        std::unique_ptr<Ort::Env> fEnv;

        /// Pointer to the ONNX inference session.
        /// This session is used to execute the model for inference.
        std::unique_ptr<Ort::Session> fSession;

        /// ONNX session options.
        /// These settings control the behavior of the inference session, such as optimization level, 
        /// execution providers, and other configuration parameters.
        Ort::SessionOptions fSessionOptions;

        /// ONNX memory info.
        /// This object provides information about memory allocation and is used during the creation of 
        /// ONNX tensors. It specifies the memory type and device (e.g., CPU, GPU).
        const OrtMemoryInfo* fInfo;
        struct MemoryInfo;

        /// Stores the input and output names for the ONNX model.
        /// These vectors contain the names of the inputs (fInames) and outputs (fOnames) that the model expects.
        std::vector<const char*> fInames;
        std::vector<const char*> fOnames;

        /// Property to specify the path to the ONNX model file.
        /// This is a configurable property that defines the location of the ONNX model file on the filesystem.
        Gaudi::Property<std::string> modelPath{this, "modelPath", "", "modelPath"};
        Gaudi::Property<double> tbeta{this, "tbeta", 0.6, "tbeta"};
        Gaudi::Property<double> td{this, "td", 0.3, "td"};

        //------------------------------------------------------------------
        //          machinery for geometry

        /// Geometry service name
        Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
        Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "uidSvc", "The name of the UniqueIDGenSvc instance"};

        /// Detector name
        Gaudi::Property<std::string> m_DCH_name{this, "DCH_name", "DCH_v2", "Name of the Drift Chamber detector"};

        /// Pointer to the geometry service
        SmartIF<IGeoSvc> m_geoSvc;

        /// Pointers to drift chamber information
        dd4hep::rec::DCH_info* dch_info;
        dd4hep::DDSegmentation::BitFieldCoder* dc_decoder;
        

};

DECLARE_COMPONENT(GGTF_tracking)