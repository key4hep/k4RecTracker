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

file(GLOB testfiles
  ${PROJECT_SOURCE_DIR}/test/*
)

file(COPY ${testfiles} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

install(TARGETS ${PackageName}
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/@{CMAKE_PROJECT_NAME}" COMPONENT dev
)

install(FILES ${testfiles} DESTINATION test)

SET(test_name "test_VTXdigitizerDetailed")
ADD_TEST(NAME ${test_name} COMMAND bash test/test_VTXdigitizerDetailed.sh)
set_tests_properties(${test_name} PROPERTIES ENVIRONMENT "K4RECTRACKER=${CMAKE_INSTALL_PREFIX}/share/k4RecTracker;PATH=${CMAKE_INSTALL_PREFIX}/bin:$ENV{PATH};PYTHONPATH=${CMAKE_INSTALL_PREFIX}/python:$ENV{PYTHONPATH};LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:${CMAKE_INSTALL_PREFIX}/lib64:$ENV{LD_LIBRARY_PATH}")
set_test_env(${test_name})
