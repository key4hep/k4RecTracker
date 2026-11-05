set(PackageName Tracking)

project(${PackageName})

include(ExternalData)
list(APPEND ExternalData_URL_TEMPLATES
  "https://key4hep.web.cern.ch:443/testFiles/k4RecTracker/%(hash)"
)

find_package(GenFit REQUIRED)
find_package(MarlinUtil REQUIRED)
list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)
FOREACH( pkg  MarlinUtil PandoraSDK )
  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()

# CONFIGURE_DEPENDS: without it a newly added .cpp is invisible to an incremental build.
file(GLOB sources CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/components/*.cpp
    ${PROJECT_SOURCE_DIR}/src/genfit_interfaces/*.cpp
    ${PROJECT_SOURCE_DIR}/src/combinatorial_track_finder/*.cpp

)

file(GLOB headers CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/include/*.h
  ${PROJECT_SOURCE_DIR}/include/genfit_interfaces/*.h
  ${PROJECT_SOURCE_DIR}/include/combinatorial_track_finder/*.h
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
  torch
  DD4hep::DDRec
  DD4hep::DDCore
  GenFit::genfit2
  Eigen3::Eigen
  ${DELPHES_LIBRARIES}
)

if(TARGET PandoraPFA::PandoraSDK)
  target_link_libraries(${PackageName} PRIVATE PandoraPFA::PandoraSDK)
endif()

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/genfit_interfaces>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/combinatorial_track_finder>
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

# Standalone unit test of the helix extrapolation maths (no Gaudi, no input files):
# analytic crossings vs an independent RK4 Lorentz-force integration.
add_executable(test_HelixExtrapolator
  ${PROJECT_SOURCE_DIR}/test/testCombinatorialTrackFinder/test_HelixExtrapolator.cpp
  ${PROJECT_SOURCE_DIR}/src/combinatorial_track_finder/HelixExtrapolator.cpp)
target_include_directories(test_HelixExtrapolator PRIVATE ${PROJECT_SOURCE_DIR}/include/combinatorial_track_finder)
target_link_libraries(test_HelixExtrapolator PRIVATE Eigen3::Eigen)
add_test(NAME test_HelixExtrapolator COMMAND test_HelixExtrapolator)

# Toy calibration of the multiple-scattering-aware kink chi2 (no Gaudi, no input files).
add_executable(test_TrackQuality
  ${PROJECT_SOURCE_DIR}/test/testCombinatorialTrackFinder/test_TrackQuality.cpp
  ${PROJECT_SOURCE_DIR}/src/combinatorial_track_finder/TrackingMath.cpp
  ${PROJECT_SOURCE_DIR}/src/combinatorial_track_finder/HelixExtrapolator.cpp)
target_include_directories(test_TrackQuality PRIVATE ${PROJECT_SOURCE_DIR}/include/combinatorial_track_finder)
target_link_libraries(test_TrackQuality PRIVATE Eigen3::Eigen)
add_test(NAME test_TrackQuality COMMAND test_TrackQuality)

# End-to-end test of CombinatorialTrackFinder: simulated muons -> digitisation -> tracks,
# then checks of the efficiency, the pT and the stored track states.
# It runs in the build directory, so its output files stay out of the source tree.
SET(test_name "test_CombinatorialTrackFinder")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testCombinatorialTrackFinder)
ADD_TEST(NAME ${test_name} COMMAND sh +x ${PROJECT_SOURCE_DIR}/test/testCombinatorialTrackFinder/test_CombinatorialTrackFinder.sh)
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/testCombinatorialTrackFinder")

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

SET(test_dir "${CMAKE_SOURCE_DIR}/Tracking/test/testTrackMerger")

# Run the setup -> k4run -> check chain once per matching mode
foreach(mode absolute significance)
  SET(test_name "test_TrackMerger_${mode}")
  SET(input_file "input_${mode}.root")
  SET(output_file "output_${mode}.root")

  ADD_TEST(NAME ${test_name}_setup
    COMMAND python3 test_track_merger.py setup --mode ${mode} --input ${input_file})
  set_test_env(${test_name}_setup)
  set_tests_properties(${test_name}_setup PROPERTIES
    WORKING_DIRECTORY "${test_dir}"
    FIXTURES_SETUP ${test_name}_input_fixture
  )

  ADD_TEST(NAME ${test_name}_run
    COMMAND k4run test_track_merger_steer.py --mode ${mode} --input ${input_file} --output ${output_file})
  set_test_env(${test_name}_run)
  set_tests_properties(${test_name}_run PROPERTIES
    WORKING_DIRECTORY "${test_dir}"
    FIXTURES_REQUIRED ${test_name}_input_fixture
    FIXTURES_SETUP ${test_name}_output_fixture
  )

  ADD_TEST(NAME ${test_name}_check
    COMMAND python3 test_track_merger.py check --mode ${mode} --output ${output_file})
  set_test_env(${test_name}_check)
  set_tests_properties(${test_name}_check PROPERTIES
    WORKING_DIRECTORY "${test_dir}"
    FIXTURES_REQUIRED ${test_name}_output_fixture
  )
endforeach()
#endif()
