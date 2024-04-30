#include "GenFitter.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <torch/torch.h>
#include "onnxruntime_cxx_api.h"

DECLARE_COMPONENT(GenFitter)

GenFitter::GenFitter(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  declareProperty("inputHits_CDC", m_input_hits_CDC, "Input CDC tracker hit collection name");
  declareProperty("inputHits_VTXIB", m_input_hits_VTXIB, "Input VTXIB tracker hit collection name");
  declareProperty("inputHits_VTXD", m_input_hits_VTXD, "Input VTXD tracker hit collection name");
  declareProperty("inputHits_VTXOB", m_input_hits_VTXOB, "Input VTXOB tracker hit collection name");
  declareProperty("outputTracks", m_output_tracks, "Output track collection name");
}

GenFitter::~GenFitter() {}

StatusCode GenFitter::initialize() { return StatusCode::SUCCESS; }

StatusCode GenFitter::execute() {
  // Get the input collection with tracker hits
  const extension::DriftChamberDigiCollection* input_hits_CDC = m_input_hits_CDC.get();
  std::cout << "Input Hit collection size CDC: " << input_hits_CDC->size() << std::endl;
  const edm4hep::TrackerHit3DCollection* inputHits_VTXD = m_input_hits_VTXD.get();
  const edm4hep::TrackerHit3DCollection* inputHits_VTXIB = m_input_hits_VTXIB.get();
  const edm4hep::TrackerHit3DCollection* inputHits_VTXOB = m_input_hits_VTXOB.get();
  std::cout << "Input Hit collection size VTXD: " << inputHits_VTXD->size() << std::endl;
  std::cout << "Input Hit collection size VTXIB: " << inputHits_VTXIB->size() << std::endl;
  std::cout << "Input Hit collection size VTXOB: " << inputHits_VTXOB->size() << std::endl;

  // // Create input list of tensors 
  // int size_CDC = input_hits_CDC->size();
  // int size_VTXD = inputHits_VTXD->size();
  // int size_VTXIB = inputHits_VTXIB->size();
  // int size_VTXOB = inputHits_VTXOB->size();

  // size_t size_total = size_CDC+size_VTXD+size_VTXIB+size_VTXOB;
  std::vector <float> ListGlobalInputs; 
 
  int it = 0;
  for (const auto& input_sim_hit : *input_hits_CDC) {
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getLeftPosition().y);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().x - input_sim_hit.getLeftPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().y - input_sim_hit.getLeftPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getRightPosition().z - input_sim_hit.getLeftPosition().z); 
    it += 1;                       
  }

  // std::cout << "Input Hit collection size CDC: " << it <<std::endl;

  for (const auto& input_sim_hit : *inputHits_VTXD) {
    ListGlobalInputs.push_back(input_sim_hit.getPosition().x);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
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
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
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
    ListGlobalInputs.push_back(input_sim_hit.getPosition().y);
    ListGlobalInputs.push_back(1.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    ListGlobalInputs.push_back(0.0);
    it += 1;                         
  }
  // std::cout << "Input Hit collection size inputHits_VTXOB: " << it <<std::endl;

  torch::Tensor input_tensor = torch::from_blob(ListGlobalInputs.data(), {it,7}, torch::kFloat32);
  // std::cout << "input_tensor: " << input_tensor <<std::endl;

  Ort::MemoryInfo memory_info = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ONNX_Runtime");
  std::cout << "Environment created" << std::endl;

  size_t total_size = it*7;
  std::vector<int64_t> tensor_shape = {it, 7};
  std::vector<Ort::Value> input_tensors;
  input_tensors.emplace_back( Ort::Value::CreateTensor<float>(memory_info, ListGlobalInputs.data(), total_size, tensor_shape.data(), tensor_shape.size()));
  // std::cout << "\ninput_tensor shape: " << print_shape(input_tensors[0].GetTensorTypeAndShapeInfo().GetShape()) << std::endl;
  // Run the model
  std::string model_path = "/afs/cern.ch/work/m/mgarciam/private/k4RecTracker_dev_0/Tracking/model_multivector_1_input.onnx";
  Ort::SessionOptions session_options;
  session_options.SetIntraOpNumThreads(1);
  session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);
  printf("Using Onnxruntime C++ API\n");
  Ort::Session session(env, model_path.c_str(), session_options);	
  printf("Starting to run inference\n");
  Ort::AllocatorWithDefaultOptions allocator;
  // get inputs and outputs
  std::vector<std::string> input_names;
  std::vector<std::int64_t> input_shapes;
  std::cout << "Input Node Name/Shape (" << input_names.size() << "):" << std::endl;
  for (std::size_t i = 0; i < session.GetInputCount(); i++) {
      input_names.emplace_back(session.GetInputName(i, allocator));
      input_shapes = session.GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
      // std::cout << "\t" << input_names.at(i) << " : " << print_shape(input_shapes) << std::endl;
  }
  for (auto& s : input_shapes) {
      if (s < 0) {
      s = 1;
      }
  }
  std::vector<std::string> output_names;
  std::cout << "Output Node Name/Shape (" << output_names.size() << "):" << std::endl;
  for (std::size_t i = 0; i < session.GetOutputCount(); i++) {
      output_names.emplace_back(session.GetOutputName(i, allocator));
      auto output_shapes = session.GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
      // std::cout << "\t" << output_names.at(i) << " : " << print_shape(output_shapes) << std::endl;
  }
  
  std::vector<const char*> input_names_char(input_names.size(), nullptr);
  std::transform(std::begin(input_names), std::end(input_names), std::begin(input_names_char),
                [&](const std::string& str) { return str.c_str(); });
                
  std::vector<const char*> output_names_char(output_names.size(), nullptr);
  std::transform(std::begin(output_names), std::end(output_names), std::begin(output_names_char),
                [&](const std::string& str) { return str.c_str(); });


  auto output_tensors = session.Run(Ort::RunOptions{nullptr}, input_names_char.data(), input_tensors.data(),
                                    input_names_char.size(), output_names_char.data(), output_names_char.size());
  std::cout << "Done!" << std::endl;
  // Get the outputs and covert to vector 
  float* floatarr = output_tensors.front().GetTensorMutableData<float>();
  std::vector<float> output_vector(floatarr, floatarr + it*4);

  torch::Tensor output_model_tensor = torch::from_blob(output_vector.data(), {it,4}, torch::kFloat32);
  std::cout << "output_model_tensor: " << output_model_tensor <<std::endl;
  // Produce dummy tracks
  edm4hep::TrackCollection* output_tracks = m_output_tracks.createAndPut();
  auto                      output_track  = output_tracks->create();
  output_track.setChi2(1.);
  output_track.setNdf(1);
  output_track.setDEdx(1.);
  return StatusCode::SUCCESS;
}

StatusCode GenFitter::finalize() { return StatusCode::SUCCESS; }
