# Header-only utilities shared between the digitisers and the tracking algorithms
add_library(k4RecTrackerUtils INTERFACE)
add_library(k4RecTracker::k4RecTrackerUtils ALIAS k4RecTrackerUtils)

target_include_directories(k4RecTrackerUtils INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(k4RecTrackerUtils INTERFACE ROOT::Hist)

install(TARGETS k4RecTrackerUtils EXPORT ${CMAKE_PROJECT_NAME}Targets)
install(DIRECTORY include/k4RecTracker DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(BUILD_TESTING)
  add_executable(testClusterCounting test/testClusterCounting.cpp)
  target_link_libraries(testClusterCounting PRIVATE k4RecTrackerUtils)
  add_test(NAME test_ClusterCounting COMMAND testClusterCounting)
endif()
