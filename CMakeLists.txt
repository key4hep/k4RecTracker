CMAKE_MINIMUM_REQUIRED(VERSION 3.12 FATAL_ERROR)

set( PackageName k4RecTracker )
project(${PackageName})

# please keep this layout for version setting, used by the automatic tagging script
set(${PackageName}_VERSION_MAJOR 0)
set(${PackageName}_VERSION_MINOR 3)
set(${PackageName}_VERSION_PATCH 0)

set(${PackageName}_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")

# dependencies
# For some reason, the find_package with DD4hep must be the first one
find_package(DD4hep REQUIRED COMPONENTS DDRec DDG4 DDParsers)
# dd4hep_set_compiler_flags() expects DD4hep_SET_RPATH to be set to ON
# otherwise it will not set the rpath when installing
set(DD4HEP_SET_RPATH ON)
dd4hep_set_compiler_flags()
find_package(ROOT COMPONENTS RIO Tree MathCore)
find_package(EDM4HEP)
find_package(k4FWCore)
find_package(Gaudi)
find_package(GSL)
#---------------------------------------------------------------

include(cmake/Key4hepConfig.cmake)
include(GNUInstallDirs)

include(CTest)
function(set_test_env _testname)
	set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT "ROOT_INCLUDE_PATH=$<$<TARGET_EXISTS:podio::podio>:$<TARGET_FILE_DIR:podio::podio>/../include>:$<$<TARGET_EXISTS:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>/../include>:$ENV{ROOT_INCLUDE_PATH}")
	set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=${PROJECT_BINARY_DIR}:${PROJECT_BINARY_DIR}/${CMAKE_PROJECT_NAME}:${PROJECT_BINARY_DIR}/genConfDir/${PackageName}:$<$<TARGET_EXISTS:ROOT::Core>:$<TARGET_FILE_DIR:ROOT::Core>>:$<$<TARGET_EXISTS:EDM4HEP::edm4hep>:$<TARGET_FILE_DIR:EDM4HEP::edm4hep>>:$<$<TARGET_EXISTS:podio::podio>:$<TARGET_FILE_DIR:podio::podio>>:$ENV{LD_LIBRARY_PATH}")
	set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}:${PROJECT_BINARY_DIR}/genConfDir:$ENV{PYTHONPATH}")
endfunction()

add_subdirectory(DCHdigi)
add_subdirectory(ARCdigi)
add_subdirectory(VTXdigi)
add_subdirectory(MUONdigi)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
add_subdirectory(Tracking)

include(cmake/CreateProjectConfig.cmake)
