#include "GenFitter.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"


torch::Tensor find_condpoints(torch::Tensor betas, torch::Tensor  unassigned, float tbeta) {
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
torch::Tensor get_clustering(std::vector<float> output_vector, int64_t num_rows,  float tbeta=0.6, float td=0.05) {

    torch::Tensor output_model_tensor = torch::from_blob(output_vector.data(), {num_rows,4}, torch::kFloat32);
    auto rows_output_model = torch::arange(0, output_model_tensor.size(0), torch::kLong);
    auto coord_ind_output_model = torch::arange(0, 3, torch::kLong);
    auto betas = output_model_tensor.index({rows_output_model,3});
    auto X = output_model_tensor.index({torch::indexing::Slice(), torch::indexing::Slice(0, 3)});
    int n_points = betas.size(0);
    auto select_condpoints = betas.gt(tbeta);
    auto indices_condpoints = find_condpoints(betas, torch::arange(n_points), tbeta);
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
        indices_condpoints = find_condpoints(betas, unassigned, tbeta);
        index_assignation += 1;
    }
    return clustering;
}

std::string print_shape(const std::vector<std::int64_t>& v) {
std::stringstream ss("");
for (std::size_t i = 0; i < v.size() - 1; i++) ss << v[i] << "x";
ss << v[v.size() - 1];
return ss.str();
}

DECLARE_COMPONENT(GenFitter)

GenFitter::GenFitter(const std::string& aName, ISvcLocator* aSvcLoc) : Gaudi::Algorithm(aName, aSvcLoc) {
  declareProperty("modelPath", modelPath, "Path to the model");
  declareProperty("inputHits_CDC", m_input_hits_CDC, "Input CDC tracker hit collection name");
  declareProperty("inputHits_VTXIB", m_input_hits_VTXIB, "Input VTXIB tracker hit collection name");
  declareProperty("inputHits_VTXD", m_input_hits_VTXD, "Input VTXD tracker hit collection name");
  declareProperty("inputHits_VTXOB", m_input_hits_VTXOB, "Input VTXOB tracker hit collection name");
  declareProperty("outputTracks", m_output_tracks, "Output track collection name");
}

GenFitter::~GenFitter() {}

StatusCode GenFitter::initialize() { 
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


StatusCode GenFitter::execute(const EventContext&) const {
  // Get the input collection with tracker hits
  const extension::DriftChamberDigiCollection* input_hits_CDC = m_input_hits_CDC.get();
  // std::cout << "Input Hit collection size CDC: " << input_hits_CDC->size() << std::endl;
  const edm4hep::TrackerHit3DCollection* inputHits_VTXD = m_input_hits_VTXD.get();
  const edm4hep::TrackerHit3DCollection* inputHits_VTXIB = m_input_hits_VTXIB.get();
  const edm4hep::TrackerHit3DCollection* inputHits_VTXOB = m_input_hits_VTXOB.get();
  // std::cout << "Input Hit collection size VTXD: " << inputHits_VTXD->size() << std::endl;
  // std::cout << "Input Hit collection size VTXIB: " << inputHits_VTXIB->size() << std::endl;
  // std::cout << "Input Hit collection size VTXOB: " << inputHits_VTXOB->size() << std::endl;

  // size_t size_total = size_CDC+size_VTXD+size_VTXIB+size_VTXOB;
  std::vector <float> ListGlobalInputs; 
 
  int it = 0;
  for (const auto& input_sim_hit : *inputHits_VTXD) {
    ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
    ListGlobalInputs.push_back(1.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0); 
    it += 1;                           
  }
  // std::cout << "Input Hit collection size inputHits_VTXD: " << it <<std::endl;
  for (const auto& input_sim_hit : *inputHits_VTXIB) {
    ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
    ListGlobalInputs.push_back(1.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0); 
    it += 1;                        
  }
  // std::cout << "Input Hit collection size inputHits_VTXIB: " << it <<std::endl;
  for (const auto& input_sim_hit : *inputHits_VTXOB) {
    ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().z);
    ListGlobalInputs.push_back(1.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    it += 1;                         
  }
  // std::cout << "Input Hit collection size inputHits_VTXOB: " << it <<std::endl;
  for (const auto& input_sim_hit : *input_hits_CDC) {
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().z);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().x - input_sim_hit.getLeftPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().y - input_sim_hit.getLeftPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().z - input_sim_hit.getLeftPosition().z); 
    it += 1;                       
  }
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
  auto clustering = get_clustering(output_vector, it);
  // std::cout << "Clustering: " << clustering << std::endl;
  // torch::Tensor output_model_tensor = torch::from_blob(output_vector.data(), {it,4}, torch::kFloat32);
  // std::cout << "output_model_tensor: " << output_model_tensor <<std::endl;
  // Produce dummy tracks
  edm4hep::TrackCollection* output_tracks = m_output_tracks.createAndPut();
  auto                      output_track  = output_tracks->create();
  output_track.setChi2(1.);
  output_track.setNdf(1);
  output_track.setDEdx(1.);
  return StatusCode::SUCCESS;
}

StatusCode GenFitter::finalize() { return StatusCode::SUCCESS; }
