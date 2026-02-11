set(PackageName VTXdigiDetailed)

project(${PackageName})

file(GLOB sources
    ${PROJECT_SOURCE_DIR}/src/*.cpp
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
  DD4hep::DDRec
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(${PackageName} PROPERTIES PUBLIC_HEADER "${headers}")

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}" COMPONENT dev
)

set(test_name "test_VTXdigitizerDetailed")
set(test_script "${CMAKE_CURRENT_SOURCE_DIR}/test/test_VTXdigitizerDetailed.sh")
add_test(NAME ${test_name}
         COMMAND bash ${test_script}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} # Le test s'exécute dans build, mais lit le script dans source
)
set_tests_properties(${test_name} PROPERTIES
    ENVIRONMENT "SOURCE_DIR_TEST=${CMAKE_CURRENT_SOURCE_DIR}/test;K4RECTRACKER=${CMAKE_INSTALL_PREFIX}/share/k4RecTracker;PATH=${CMAKE_INSTALL_PREFIX}/bin:$ENV{PATH};PYTHONPATH=${CMAKE_INSTALL_PREFIX}/python:${CMAKE_CURRENT_SOURCE_DIR}/test:$ENV{PYTHONPATH};LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib64:$ENV{LD_LIBRARY_PATH}"
)
if(COMMAND set_test_env)
  set_test_env(${test_name})
endif()