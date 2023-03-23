CMAKE_MINIMUM_REQUIRED(VERSION 3.12)


project(k4ProjectTemplate)

# please keep this layout for version setting, used by the automatic tagging script
set(PACKAGE_VERSION_MAJOR 1)
set(PACKAGE_VERSION_MINOR 2)
set(PACKAGE_VERSION_PATCH 3)

set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

# dependencies are declared like this
find_package(ROOT COMPONENTS RIO Tree)
find_package(EDM4HEP)
find_package(k4FWCore)
find_package(Gaudi)
#---------------------------------------------------------------

include(GNUInstallDirs)

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

if(NOT CMAKE_CXX_STANDARD MATCHES "14|17|20")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()

include(CTest)

add_subdirectory(DCHdigi)

include(cmake/CreateProjectConfig.cmake)
