###############################################################################
# This file contains the boilerplate to create the package config.cmake file
# used by other packages to find this package
# In principle no modifications are needed, except for the dependencies used
# in which case please modify the Project_Config_Template.cmake.in file
# and in rare cases the parts related to PATH_VARS here and in
# Project_Config_Template.cmake.in, but this should be rarely needed, because
# this information should go via the TARGETS
#
# This setup requires the use of GNUInstallDirs
###############################################################################

configure_file(${CMAKE_CURRENT_LIST_DIR}/package_version_template.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/include/${CMAKE_PROJECT_NAME}/${CMAKE_PROJECT_NAME}Version.h)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  cmake/Project_Config_Template.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}
  # Adapt if other variables are needed, also modify Project_Config_template.cmake.in
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}
  COMPATIBILITY SameMajorVersion )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME} )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/${CMAKE_PROJECT_NAME}/${CMAKE_PROJECT_NAME}Version.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME} )

install(EXPORT ${CMAKE_PROJECT_NAME}Targets
  NAMESPACE ${CMAKE_PROJECT_NAME}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}"
  )
