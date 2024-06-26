#include "GenFitter_evaluation.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"
#include <ATen/ATen.h>
#include "edm4hep/TrackerHit3D.h"
#include "extension/MutableTrackerHit3D.h"

torch::Tensor find_condpoints1(torch::Tensor betas, torch::Tensor  unassigned, float tbeta) {
    int n_points = unassigned.size(0);
    int size_b = betas.size(0);
    auto select_condpoints = betas.gt(tbeta);
    auto mask_unassigned = torch::zeros({size_b}, torch::dtype(torch::kBool));
    for (int i = 0; i < n_points-1; ++i) {
        auto ii = unassigned.index({i});
        mask_unassigned.index_put_({ii.squeeze()}, true);
    }
    select_condpoints = mask_unassigned * select_condpoints;
    auto indices_condpoints = select_condpoints.nonzero();
    auto betas_condpoints = -betas.index({indices_condpoints});
    auto sorted_indices = torch::argsort(betas_condpoints, /*dim=*/0, /*descending=*/false);
    indices_condpoints = indices_condpoints.index({sorted_indices});
    return indices_condpoints;
}

// Main clustering function
torch::Tensor get_clustering1(std::vector<float> output_vector, int64_t num_rows,  float tbeta=0.6, float td=0.05) {

    torch::Tensor output_model_tensor = torch::from_blob(output_vector.data(), {num_rows,4}, torch::kFloat32);
    auto rows_output_model = torch::arange(0, output_model_tensor.size(0), torch::kLong);
    auto coord_ind_output_model = torch::arange(0, 3, torch::kLong);
    auto betas = output_model_tensor.index({rows_output_model,3});
    auto X = output_model_tensor.index({torch::indexing::Slice(), torch::indexing::Slice(0, 3)});
    int n_points = betas.size(0);
    auto select_condpoints = betas.gt(tbeta);
    auto indices_condpoints = find_condpoints1(betas, torch::arange(n_points), tbeta);
    auto clustering = -1 * torch::ones({n_points}, torch::kLong);
    int index_assignation = 0;
    auto unassigned = torch::arange(n_points);
    while (indices_condpoints.size(0) > 0 && unassigned.size(0) > 0) {
        auto index_condpoint = indices_condpoints.index({0});
        auto d = (X.index({unassigned}) - X.index({index_condpoint})).norm(2, -1).squeeze();
        auto mask_distance = d.lt(td);
        auto  mask_distance_ind = mask_distance.nonzero();
        auto assigned_to_this_condpoint = unassigned.index({mask_distance_ind});
        clustering.index_put_({assigned_to_this_condpoint}, index_assignation);
        auto mask_distance_out = d.ge(td);
        auto  mask_distance_ind_out = mask_distance_out.squeeze(0).nonzero();
        unassigned = unassigned.index({mask_distance_ind_out});
        indices_condpoints = find_condpoints1(betas, unassigned, tbeta);
        index_assignation += 1;
    }
    return clustering;
}

std::string print_shape1(const std::vector<std::int64_t>& v) {
std::stringstream ss("");
for (std::size_t i = 0; i < v.size() - 1; i++) ss << v[i] << "x";
ss << v[v.size() - 1];
return ss.str();
}

DECLARE_COMPONENT(GenFitter_eval)

GenFitter_eval::GenFitter_eval(const std::string& aName, ISvcLocator* aSvcLoc) : Gaudi::Algorithm(aName, aSvcLoc) {
  declareProperty("modelPath", modelPath, "Path to the model");
  declareProperty("inputHits_CDC", m_input_hits_CDC, "Input CDC tracker hit collection name");
  declareProperty("inputHits_VTXIB", m_input_hits_VTXIB, "Input VTXIB tracker hit collection name");
  declareProperty("inputHits_VTXD", m_input_hits_VTXD, "Input VTXD tracker hit collection name");
  declareProperty("inputHits_VTXOB", m_input_hits_VTXOB, "Input VTXOB tracker hit collection name");
  declareProperty("inputAssociation_CDC_sim", m_input_Association_CDC, "Input CDC association  collection name");
  declareProperty("inputHits_CDC_sim", m_input_hits_CDC_sim, "Input CDC sim hit collection name");
  declareProperty("inputHits_VTXIB_sim", m_input_hits_VTXIB_sim, "Input VTXIB sim hit collection name");
  declareProperty("inputHits_VTXD_sim", m_input_hits_VTXD_sim, "Input VTXD sim hit collection name");
  declareProperty("inputHits_VTXOB_sim", m_input_hits_VTXOB_sim, "Input VTXOB sim hit collection name");
  declareProperty("outputTracks", m_output_tracks, "Output track collection name");
  declareProperty("outputHits", m_output_hits, "Output hits collection name");
}

GenFitter_eval::~GenFitter_eval() {}

StatusCode GenFitter_eval::initialize() { 
  fInfo = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  auto envLocal = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
  fEnv          = std::move(envLocal);
  // std::cout << "Environment created" << std::endl;
  // std::string model_path = "/afs/cern.ch/work/m/mgarciam/private/k4RecTracker_dev_0/Tracking/model_multivector_1_input.onnx";
  fSessionOptions.SetIntraOpNumThreads(1);
  fSessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);
  // printf("Using Onnxruntime C++ API\n");
  auto sessionLocal = std::make_unique<Ort::Session>(*fEnv, modelPath.c_str(), fSessionOptions);	
  fSession          = std::move(sessionLocal);
  // printf("Starting to run inference\n");
  Ort::AllocatorWithDefaultOptions allocator;
  std::size_t i = 0;
  const auto input_name = fSession->GetInputNameAllocated(i, allocator).release();
  
  const auto output_names = fSession->GetOutputNameAllocated(i, allocator).release();
  
  fInames.push_back(input_name);
  fOnames.push_back(output_names);
  return StatusCode::SUCCESS;}


StatusCode GenFitter_eval::execute(const EventContext&) const {
  // Get the input collection with tracker hits
  const extension::DriftChamberDigiCollection* input_hits_CDC = m_input_hits_CDC.get();
  // std::cout << "Input Hit collection size CDC: " << input_hits_CDC->size() << std::endl;
  const edm4hep::TrackerHit3DCollection* inputHits_VTXD = m_input_hits_VTXD.get();
  const edm4hep::TrackerHit3DCollection* inputHits_VTXIB = m_input_hits_VTXIB.get();
  const edm4hep::TrackerHit3DCollection* inputHits_VTXOB = m_input_hits_VTXOB.get();

  const edm4hep::SimTrackerHitCollection* inputHits_VTXD_sim = m_input_hits_VTXD_sim.get();
  const edm4hep::SimTrackerHitCollection* inputHits_VTXIB_sim = m_input_hits_VTXOB_sim.get();
  const edm4hep::SimTrackerHitCollection* inputHits_VTXOB_sim = m_input_hits_VTXOB_sim.get();
  const edm4hep::SimTrackerHitCollection* inputHits_CDC_sim = m_input_hits_CDC_sim.get();
  const extension::MCRecoDriftChamberDigiAssociation* inputAssociation_CDC_sim = m_input_Association_CDC.get();
  std::cout << "Input Hit collection size VTXD: " << inputHits_VTXD_sim->size() << std::endl;
  // std::cout << "Input Hit collection size VTXIB: " << inputHits_VTXIB->size() << std::endl;
  // std::cout << "Input Hit collection size VTXOB: " << inputHits_VTXOB->size() << std::endl;
  
  // For now there is a 1to1 correspondence between sim and digi so we can get the list of mc for each hit
  std::vector <float> ListHitMC_VTXD; 
  for (const auto& input_sim_hit : *inputHits_VTXD_sim) {
    auto MC_particle = input_sim_hit.getParticle();
    auto object_id_MC = MC_particle.getObjectID();
    auto index_MC = object_id_MC.index;
    ListHitMC_VTXD.push_back(index_MC);                 
  }
  std::vector <float> ListHitMC_VTXIB; 
  for (const auto& input_sim_hit : *inputHits_VTXIB_sim) {
    auto MC_particle = input_sim_hit.getParticle();
    auto object_id_MC = MC_particle.getObjectID();
    auto index_MC = object_id_MC.index;
    ListHitMC_VTXIB.push_back(index_MC);                 
  }
  std::vector <float> ListHitMC_VTXOB; 
  for (const auto& input_sim_hit : *inputHits_VTXOB_sim) {
    auto MC_particle = input_sim_hit.getParticle();
    auto object_id_MC = MC_particle.getObjectID();
    auto index_MC = object_id_MC.index;
    ListHitMC_VTXOB.push_back(index_MC);                 
  }

  // size_t size_total = size_CDC+size_VTXD+size_VTXIB+size_VTXOB;
  std::vector <float> ListGlobalInputs; 
  std::vector <float> ListHitType_VTXD;
  int it = 0;
  int it_0 = 0;
  for (const auto& input_sim_hit : *inputHits_VTXD) {
    ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
    ListGlobalInputs.push_back(1.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0); 
    ListHitType_VTXD.push_back(it);
    it += 1;  
    it_0 += 1;                        
  }
  torch::Tensor ListHitType_VTXD_tensor = torch::from_blob(ListHitType_VTXD.data(), {it_0}, torch::kFloat32);
  std::cout << "Here1" << std::endl;
  std::vector <float> ListHitType_VTXIB;
  int it_1 = 0;
  for (const auto& input_sim_hit : *inputHits_VTXIB) {
    ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
    ListGlobalInputs.push_back(1.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0); 
    ListHitType_VTXIB.push_back(it);
    it += 1; 
    it_1 += 1;                      
  }
  torch::Tensor ListHitType_VTXIB_tensor = torch::from_blob(ListHitType_VTXIB.data(), {it_1}, torch::kFloat32);
  // std::cout << "Input Hit collection size inputHits_VTXIB: " << it <<std::endl;
  std::vector <float> ListHitType_VTXOB;
  int it_2 = 0;
  for (const auto& input_sim_hit : *inputHits_VTXOB) {
    ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
    ListGlobalInputs.push_back(1.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    ListHitType_VTXOB.push_back(it);
    it += 1;  
    it_2 += 1;                   
  }
  torch::Tensor ListHitType_VTXOB_tensor = torch::from_blob(ListHitType_VTXOB.data(), {it_2}, torch::kFloat32);
  // std::cout << "Input Hit collection size inputHits_VTXOB: " << it <<std::endl;
  int it_3 = 0;

  std::vector <float> ListHitMC_CDC; 
  std::vector <float> ListHitType_CDC;
  for (const auto& input_sim_hit : *input_hits_CDC) {
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().z);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().x - input_sim_hit.getLeftPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().y - input_sim_hit.getLeftPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().z - input_sim_hit.getLeftPosition().z); 
    ListHitType_CDC.push_back(it);

    // auto association_cdc = inputAssociation_CDC_sim[it_3];
    // auto hit_reference = association_cdc.getSim();
    // auto corresponding_sim_hit_CDC = hit_reference;
    // auto corresponding_sim_hit_CDC_id = corresponding_sim_hit_CDC.getObjectID();
    // std::cout << "Input Hit association CDC: " << corresponding_sim_hit_CDC_id.index <<std::endl;
    // auto corresponding_sim_hit_CDC = inputHits_CDC_sim->at(index_hit);
    // auto MC_particle = corresponding_sim_hit_CDC.getParticle();
    // auto object_id_MC = MC_particle.getObjectID();
    // auto index_MC = object_id_MC.index;
    // ListHitMC_CDC.push_back(index_MC);
    it += 1;    
    it_3 += 1;                     
  }
  torch::Tensor ListHitType_CDC_tensor = torch::from_blob(ListHitType_CDC.data(), {it_3}, torch::kFloat32);
  // std::cout << "Input Hit collection size CDC: " << it <<std::endl;
  // torch::Tensor input_tensor = torch::from_blob(ListGlobalInputs.data(), {it,7}, torch::kFloat32);
  // std::cout << "input_tensor: " << input_tensor <<std::endl;
  size_t total_size = it*7;
  std::vector<int64_t> tensor_shape = {it, 7};
  std::vector<Ort::Value> input_tensors;
  input_tensors.emplace_back( Ort::Value::CreateTensor<float>(fInfo, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));

  auto output_tensors = fSession->Run(Ort::RunOptions{nullptr}, fInames.data(), input_tensors.data(),
                                    fInames.size(), fOnames.data(), fOnames.size());
  // std::cout << "Done!" << std::endl;
  // Get the outputs and covert to vector 
  float* floatarr = output_tensors.front().GetTensorMutableData<float>();
  std::vector<float> output_vector(floatarr, floatarr + it*4);
  auto clustering = get_clustering1(output_vector, it);
  torch::Tensor unique_tensor;
  torch::Tensor inverse_indices;
  std::tie(unique_tensor, inverse_indices)  = at::_unique(clustering, true, true);
  
  // std::cout << "Clustering: " << clustering << std::endl;
  // std::cout << "inverse_indices: " << inverse_indices << std::endl;

  
  // std::cout << "indices: " << clustering.index({indices}) << std::endl;
  extension::TrackCollection* output_tracks = m_output_tracks.createAndPut();
  extension::TrackerHit3DCollection* output_hits= m_output_hits.createAndPut();
  int64_t number_of_tracks = unique_tensor.numel(); 
  std::cout << "number_of_tracks: " << number_of_tracks << std::endl;
  std::cout << "Writing data" << std::endl;
  for (int i = 0; i < number_of_tracks; ++i) {
      auto id_of_track = unique_tensor.index({i});
      auto output_track  = output_tracks->create();
      // set global properties
      output_track.setChi2(1.);
      output_track.setNdf(1);
      output_track.setDEdx(1.);
      torch::Tensor mask = (clustering == id_of_track);
      torch::Tensor indices = torch::nonzero(mask);
      int64_t number_of_hits = indices.numel();
      for (int j = 0; j < number_of_hits; ++j) {
        auto index_id = indices.index({j});
        torch::Tensor mask_VTXD = (ListHitType_VTXD_tensor == index_id);
        torch::Tensor mask_VTXIB = (ListHitType_VTXIB_tensor == index_id);
        torch::Tensor mask_VTOB = (ListHitType_VTXOB_tensor == index_id);
        torch::Tensor mask_CDC = (ListHitType_CDC_tensor == index_id);
        if ((torch::sum(mask_VTXD)>0).item<bool>()){
          // The hit belong to vtxd
          auto hit = inputHits_VTXD->at(index_id.item<int>());
          auto hit_extension  = output_hits->create();
          hit_extension.setCellID(hit.getCellID());
          hit_extension.setType(1);
          hit_extension.setEDep(ListHitMC_VTXD[index_id.item<int>()]);
          hit_extension.setPosition(hit.getPosition());
          // output_track.addToTrackerHits(hit_extension);
          
          // add this hit to some collection
        } else if ((torch::sum(mask_VTXIB)>0).item<bool>()){
          index_id = index_id-it_0;
          auto hit = inputHits_VTXIB->at(index_id.item<int>());
          auto hit_extension  = output_hits->create();
          hit_extension.setCellID(hit.getCellID());
          hit_extension.setType(1);
          // hit_extension.setEDep(ListHitMC_VTXIB[index_id.item<int>()]);
          hit_extension.setPosition(hit.getPosition());
          // output_track.addToTrackerHits(hit_extension);
        } else if ((torch::sum(mask_VTOB)>0).item<bool>()){
          index_id = index_id-(it_1+it_0);
          auto hit = inputHits_VTXOB->at(index_id.item<int>());
          auto hit_extension  = output_hits->create();
          hit_extension.setCellID(hit.getCellID());
          hit_extension.setType(1);
          // hit_extension.setEDep(ListHitMC_VTXOB[index_id.item<int>()]);
          hit_extension.setPosition(hit.getPosition());
          output_track.addToTrackerHits(hit_extension);
        } else if ((torch::sum(mask_CDC)>0).item<bool>()){
          index_id = index_id-(it_1+it_2 +it_0);
          auto hit = input_hits_CDC->at(index_id.item<int>());
          auto hit_extension  = output_hits->create();
          hit_extension.setCellID(hit.getCellID());
          hit_extension.setType(0);
          //hit_extension.setEDep(ListHitMC_CDC[index_id.item<int>()]);
          // hit_extension.setEDep(0);
          hit_extension.setPosition(hit.getLeftPosition());
          output_track.addToTrackerHits(hit_extension);
        }
    }
  }
  
 
  return StatusCode::SUCCESS;
}

StatusCode GenFitter_eval::finalize() { return StatusCode::SUCCESS; }

