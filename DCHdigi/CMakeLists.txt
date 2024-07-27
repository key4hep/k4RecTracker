set(PackageName DCHdigi)

project(${PackageName})

# Build the extension data model and link it against the upstream model
PODIO_GENERATE_DATAMODEL(extension dataFormatExtension/driftChamberHit.yaml ext_headers ext_sources
  UPSTREAM_EDM edm4hep:${EDM4HEP_DATA_DIR}/edm4hep.yaml
  IO_BACKEND_HANDLERS ${PODIO_IO_HANDLERS}
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
PODIO_ADD_DATAMODEL_CORE_LIB(extension "${ext_headers}" "${ext_sources}"
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(extension PUBLIC EDM4HEP::edm4hep)
PODIO_ADD_ROOT_IO_DICT(extensionDict extension "${ext_headers}" src/selection.xml
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
#PODIO_ADD_SIO_IO_BLOCKS(extension "${ext_headers}" "${ext_sources}")
add_library(extension::extensionDict ALIAS extensionDict)
list(APPEND EXTENSION_INSTALL_LIBS extension extensionDict)
install(TARGETS ${EXTENSION_INSTALL_LIBS}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/extension"
  COMPONENT dev)

install(FILES
  "${PROJECT_BINARY_DIR}/extensionDictDict.rootmap"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)

install(FILES
  dataFormatExtension/driftChamberHit.yaml
  DESTINATION "${CMAKE_INSTALL_DATADIR}/extension" COMPONENT dev)

install(FILES
  "${PROJECT_BINARY_DIR}/libextensionDict_rdict.pcm"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)

file(GLOB sources
    ${PROJECT_SOURCE_DIR}/src/*.cpp
)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

gaudi_add_module(${PackageName}
  SOURCES ${sources}
  LINK
  Gaudi::GaudiAlgLib
  Gaudi::GaudiKernel
  EDM4HEP::edm4hep
  extensionDict
  k4FWCore::k4FWCore
  k4FWCore::k4Interface
  DD4hep::DDRec
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${CMAKE_INSTALL_PREFIX}/extension
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
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/@{CMAKE_PROJECT_NAME}" COMPONENT dev
)

install(FILES ${scripts} DESTINATION test)

SET(test_name "test_DCHsimpleDigitizer")
ADD_TEST(NAME ${test_name} COMMAND k4run test/runDCHsimpleDigitizer.py)
set_test_env(${test_name})

SET(test_name "test_DCHsimpleDigitizerExtendedEdm")
ADD_TEST(NAME ${test_name} COMMAND k4run test/runDCHsimpleDigitizerExtendedEdm.py)
set_test_env(${test_name})
