CMAKE_MINIMUM_REQUIRED(VERSION 3.12 FATAL_ERROR)

set( PackageName k4RecTracker )
project(${PackageName} LANGUAGES CXX)

# please keep this layout for version setting, used by the automatic tagging script
set(${PackageName}_VERSION_MAJOR 0)
set(${PackageName}_VERSION_MINOR 5)
set(${PackageName}_VERSION_PATCH 0)

set(${PackageName}_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

# Let CMake know where to find FindXXX.cmake files
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# dependencies
# For some reason, the find_package with DD4hep must be the first one
find_package(DD4hep REQUIRED COMPONENTS DDRec DDG4 DDParsers)
# dd4hep_set_compiler_flags() expects DD4hep_SET_RPATH to be set to ON
# otherwise it will not set the rpath when installing
set(DD4HEP_SET_RPATH ON)
dd4hep_set_compiler_flags()
find_package(ROOT COMPONENTS RIO Tree MathCore)
find_package(podio 1.2 REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(k4FWCore 1.3 REQUIRED)
find_package(Gaudi REQUIRED)
find_package(Delphes REQUIRED)
find_package(GSL REQUIRED)
find_package(PandoraSDK REQUIRED)
#---------------------------------------------------------------

find_package(onnxruntime REQUIRED)
if(onnxruntime_FOUND)
  message("ONNX Runtime inference library found.")
else()
  message("ONNX Runtime not found!")
endif()


find_package(Torch REQUIRED)
if(Torch_FOUND)
  message("Torch inference library found.")
else()
  message("Torch not found!")
endif()

#---------------------------------------------------------------

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 20 CACHE STRING "")

if(NOT CMAKE_CXX_STANDARD MATCHES "14|17|20")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

#---------------------------------------------------------------

include(cmake/Key4hepConfig.cmake)
include(GNUInstallDirs)

include(CTest)
function(set_test_env _testname)
  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
    "ROOT_INCLUDE_PATH=$<$<TARGET_EXISTS:podio::podio>:$<TARGET_FILE_DIR:podio::podio>/../include>:$<$<TARGET_EXISTS:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>/../include>:$ENV{ROOT_INCLUDE_PATH}"
  )

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
    "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/.plugins:${CMAKE_BINARY_DIR}/DCHdigi:${CMAKE_INSTALL_PREFIX}/lib64:${PROJECT_BINARY_DIR}:${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}:${PROJECT_BINARY_DIR}/genConfDir/${PackageName}:${CMAKE_BINARY_DIR}/DCHdigi/genConfDir/DCHdigi:$<$<TARGET_EXISTS:ROOT::Core>:$<TARGET_FILE_DIR:ROOT::Core>>:$<$<TARGET_EXISTS:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>>:$<$<TARGET_EXISTS:podio::podio>:$<TARGET_FILE_DIR:podio::podio>>:$ENV{LD_LIBRARY_PATH}"
  )

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
    "PYTHONPATH=${CMAKE_BINARY_DIR}:${PROJECT_BINARY_DIR}/genConfDir:${CMAKE_BINARY_DIR}/DCHdigi/genConfDir:$ENV{PYTHONPATH}"
  )

  set_property(TEST ${_testname} APPEND PROPERTY ENVIRONMENT
    "GAUDI_PLUGIN_PATH=${CMAKE_BINARY_DIR}/.plugins:${CMAKE_BINARY_DIR}/DCHdigi:${CMAKE_INSTALL_PREFIX}/lib64:$ENV{GAUDI_PLUGIN_PATH}"
  )
endfunction()


add_subdirectory(DCHdigi)
add_subdirectory(ARCdigi)
add_subdirectory(VTXdigi)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
add_subdirectory(Tracking)

include(cmake/CreateProjectConfig.cmake)