set(PackageName Tracking)

project(${PackageName} LANGUAGES CXX)

find_package(GenFit REQUIRED)
find_package(MarlinUtil REQUIRED)
find_package(k4Reco REQUIRED)
list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)
FOREACH( pkg  MarlinUtil PandoraSDK k4Reco )
  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()

file(GLOB sources
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/components/*.cpp
    ${PROJECT_SOURCE_DIR}/src/genfit_interfaces/*.cpp

)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
  ${PROJECT_SOURCE_DIR}/include/genfit_interfaces/*.h
)

gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
  Gaudi::GaudiKernel
  EDM4HEP::edm4hep
  k4FWCore::k4FWCore
  k4FWCore::k4Interface
  k4geo::detectorCommon
  extensionDict
  onnxruntime::onnxruntime
  DD4hep::DDRec
  DD4hep::DDCore
  GenFit::genfit2
  Eigen3::Eigen
  ${DELPHES_LIBRARIES}
  k4Reco::GaudiTrkUtils
)

if(TARGET PandoraPFA::PandoraSDK)
  target_link_libraries(${PackageName} PRIVATE PandoraPFA::PandoraSDK)
endif()

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/genfit_interfaces>
  $<INSTALL_INTERFACE:include>
  ${MarlinUtil_INCLUDE_DIRS}
  ${DELPHES_INCLUDE_DIRS}
)

set_target_properties(${PackageName} PROPERTIES PUBLIC_HEADER "${headers}")

file(GLOB scripts
  ${PROJECT_SOURCE_DIR}/test/*.py
)

file(COPY ${scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}" COMPONENT dev
)

install(FILES ${scripts} DESTINATION test)

if(BUILD_TESTING)
  include(ExternalData)
  list(APPEND ExternalData_URL_TEMPLATES
    "https://key4hep.web.cern.ch:443/testFiles/k4RecTracker/%(hash)"
  )

  SET(test_name "test_TracksFromGenParticles")
  ADD_TEST(NAME ${test_name} COMMAND sh +x test_TracksFromGenParticles.sh)
  set_test_env(${test_name})
  set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/Tracking/test/testTracksFromGenParticles")

  SET(test_name "test_TrackFinder")
  ExternalData_Add_Test(test_TrackFinder
          NAME ${test_name}
          COMMAND sh +x Tracking/test/testTrackFinder/test_trackFinder.sh DATA{${CMAKE_SOURCE_DIR}/Tracking/test/inputFiles/SimpleGatrIDEAv3o1.onnx})
  set_test_env(${test_name})
  set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
  ExternalData_Add_Target(${test_name})

  SET(test_name "test_TrackFitter")
  ADD_TEST(NAME ${test_name} COMMAND sh +x test_trackFitter.sh)
  set_test_env(${test_name})
  set_tests_properties(${test_name} PROPERTIES DEPENDS test_TrackFinder)
  set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/Tracking/test/testTrackFitter")

  set(ILDConfig_DIR "${CMAKE_CURRENT_BINARY_DIR}/ILDConfig/StandardConfig/production")
  ADD_TEST(NAME "clone_ILDConfig"
    COMMAND bash -c "if ! test -d ILDConfig; then git clone https://github.com/iLCSoft/ILDConfig.git --depth 1; fi"
  )
  set_tests_properties(clone_ILDConfig
    PROPERTIES
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
      FIXTURES_SETUP ensure_ILDConfig
  )

  # TrackMerger tests: simulate gun muons for each detector model, then run the
  # ILD reconstruction (tracking only) on them and verify the merged tracks
  foreach(det_model IN ITEMS ILD_FCCee_v01 ILD_FCCee_v02)
    SET(test_name "test_TrackMerger_sim_${det_model}")
    ADD_TEST(NAME ${test_name}
      COMMAND ddsim --enableGun --gun.distribution "cos(theta)" --gun.energy "5*GeV" --gun.particle "mu-" --gun.multiplicity 10
                    --numberOfEvents 10 --random.enableEventSeed --random.seed 42
                    --outputFile ${det_model}_SIM.edm4hep.root
                    --steeringFile ddsim_steer.py
                    --compactFile $ENV{K4GEO}/FCCee/ILD_FCCee/compact/${det_model}/${det_model}.xml
    )
    set_test_env(${test_name})
    set_tests_properties(${test_name}
      PROPERTIES
        WORKING_DIRECTORY "${ILDConfig_DIR}"
        FIXTURES_REQUIRED ensure_ILDConfig
        FIXTURES_SETUP TrackMerger_sim_${det_model}
    )

    SET(test_name "test_TrackMerger_reco_${det_model}")
    ADD_TEST(NAME ${test_name}
      COMMAND bash +x ${CMAKE_CURRENT_SOURCE_DIR}/test/testTrackMerger/test_track_merger.sh ${det_model}
    )
    set_test_env(${test_name})
    set_tests_properties(${test_name}
      PROPERTIES
        WORKING_DIRECTORY "${ILDConfig_DIR}"
        FIXTURES_REQUIRED "ensure_ILDConfig;TrackMerger_sim_${det_model}"
    )
  endforeach()
endif()
