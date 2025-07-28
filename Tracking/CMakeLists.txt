set(PackageName Tracking)

project(${PackageName})

include(ExternalData)
list(APPEND ExternalData_URL_TEMPLATES
  "https://key4hep.web.cern.ch:443/testFiles/k4RecTracker/%(hash)"
)

#find_package(GenFit)
find_package(MarlinUtil REQUIRED)
list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)
FOREACH( pkg  MarlinUtil PandoraSDK )
  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()
#if (GenFit_FOUND)

file(GLOB sources
    ${PROJECT_SOURCE_DIR}/src/*.cpp
    ${PROJECT_SOURCE_DIR}/components/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
  Gaudi::GaudiKernel
  EDM4HEP::edm4hep
  k4FWCore::k4FWCore
  k4FWCore::k4Interface
  extensionDict
  onnxruntime::onnxruntime
  torch
  DD4hep::DDRec
  DD4hep::DDCore
  #GenFit::genfit2
  ${DELPHES_LIBRARIES}
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
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

SET(test_name "test_TracksFromGenParticles")
ADD_TEST(NAME ${test_name} COMMAND sh +x test_TracksFromGenParticles.sh)
set_test_env(${test_name})
set_tests_properties(${test_name} PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/Tracking/test/testTracksFromGenParticles")

ExternalData_Add_Test(test_TrackFinder
        NAME test_TrackFinder
        COMMAND COMMAND sh +x Tracking/test/testTrackFinder/test_trackFinder.sh DATA{${CMAKE_SOURCE_DIR}/Tracking/test/inputFiles/SimpleGatrIDEAv3o1.onnx})
set_test_env(test_TrackFinder)
set_tests_properties( test_TrackFinder PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
ExternalData_Add_Target(test_TrackFinder)

#endif()
