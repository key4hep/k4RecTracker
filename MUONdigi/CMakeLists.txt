set(PackageName MUONdigi)

project(${PackageName})

# Build the extension data model and link it against the upstream model
PODIO_GENERATE_DATAMODEL(muon_extension dataFormatExtension/muonSystemHit.yaml ext_headers ext_sources
  UPSTREAM_EDM edm4hep:${EDM4HEP_DATA_DIR}/edm4hep.yaml
  IO_BACKEND_HANDLERS ${PODIO_IO_HANDLERS}
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
PODIO_ADD_DATAMODEL_CORE_LIB(muon_extension "${ext_headers}" "${ext_sources}"
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(muon_extension PUBLIC EDM4HEP::edm4hep)
PODIO_ADD_ROOT_IO_DICT(muon_extensionDict muon_extension "${ext_headers}" src/selection.xml
  OUTPUT_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
#PODIO_ADD_SIO_IO_BLOCKS(muon_extension "${ext_headers}" "${ext_sources}")
add_library(muon_extension::muon_extensionDict ALIAS muon_extensionDict)
list(APPEND EXTENSION_INSTALL_LIBS muon_extension muon_extensionDict)
install(TARGETS ${EXTENSION_INSTALL_LIBS}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/muon_extension"
  COMPONENT dev)

install(FILES
  "${PROJECT_BINARY_DIR}/muon_extensionDictDict.rootmap"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT dev)

install(FILES
  dataFormatExtension/muonSystemHit.yaml
  DESTINATION "${CMAKE_INSTALL_DATADIR}/muon_extension" COMPONENT dev)

install(FILES
  "${PROJECT_BINARY_DIR}/libmuon_extensionDict_rdict.pcm"
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
  muon_extensionDict
  k4FWCore::k4FWCore
  DD4hep::DDRec
)

target_include_directories(${PackageName} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${CMAKE_INSTALL_PREFIX}/muon_extension
)

set_target_properties(MUONdigi PROPERTIES PUBLIC_HEADER "${headers}")

file(GLOB scripts
  ${PROJECT_SOURCE_DIR}/test/*.py
)

file(COPY ${scripts} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test)

install(TARGETS MUONdigi
  EXPORT ${CMAKE_PROJECT_NAME}Targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/@{CMAKE_PROJECT_NAME}" COMPONENT dev
)

install(FILES ${scripts} DESTINATION test)

SET(test_name "test_MUONsimpleDigitizer")
ADD_TEST(NAME t_${test_name} COMMAND k4run test/runMUONsimpleDigitizer.py)