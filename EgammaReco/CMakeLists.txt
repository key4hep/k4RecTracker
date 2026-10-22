set(PackageName EgammaReco)

project(${PackageName})

find_package(MarlinUtil REQUIRED)
list(APPEND CMAKE_MODULE_PATH $ENV{PANDORAPFA}/cmakemodules)
FOREACH( pkg  MarlinUtil PandoraSDK )
  IF( ${pkg}_FOUND )
    INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
    ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()

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
                 k4FWCore::k4FWCore
                 EDM4HEP::edm4hep
                 EDM4HEP::edm4hep
                 k4FWCore::k4Interface
                 k4geo::detectorCommon
                 DD4hep::DDRec
                 DD4hep::DDCore
                 ${MarlinUtil_LIBRARIES}
                 ${PandoraSDK_LIBRARIES}
                 )

target_include_directories(${PackageName} PUBLIC
                   ${MarlinUtil_INCLUDE_DIRS}
                   ${PandoraSDK_INCLUDE_DIRS}
                   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                   $<INSTALL_INTERFACE:include>
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
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}"COMPONENT dev)

install(FILES ${scripts} DESTINATION test)

SET(test_name "test_ElectronDecorator")
ADD_TEST(NAME ${test_name} COMMAND k4run test/ElectronDecorator.py)
set_test_env(${test_name})
