# The name of our project is "COSM". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${COSM_SOURCE_DIR} and to the root binary directory of the
# project as ${COSM_BINARY_DIR}.
#
PROJECT(COSM)
cmake_minimum_required(VERSION 2.4)
INCLUDE(InstallRequiredSystemLibraries)

if(APPLE)
   set(__APPLE_QUARTZ__ 1)
   set(HAVE_STRTOLL 1)
   set(HAVE_STRCASECMP 1)
   set(HAVE_DIRENT_H 1)
   set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -framework AGL -framework Cocoa")
endif(APPLE)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Computational Optical Sectioning Microscopy (COSM)")
SET(CPACK_PACKAGE_VENDOR "Prezise Solutions")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "9")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "cosm-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
#SET(CPACK_INSTALL_CMAKE_PROJECTS "/Users/einirvaldimarsson/xcosm/build/debug/;COSM;ALL;/")

IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
#  SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
#  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\MyExecutable.exe")
#  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} COSM")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.prezise.com";)
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\www.prezise.com";)
  SET(CPACK_NSIS_CONTACT "support@prezise.com")
  SET(CPACK_NSIS_MODIFY_PATH ON)
ELSE(WIN32 AND NOT UNIX)
#  SET(CPACK_STRIP_FILES "/bin/Viewer;/bin/PsfGeneratorGUI;/bin/EstimationGui;/bin/ToolsGui")
  SET(CPACK_SOURCE_STRIP_FILES "")
ENDIF(WIN32 AND NOT UNIX)

SET(CPACK_PACKAGE_EXECUTABLES "Viewer;COSM 3D Viewer" "PsfGeneratorGui;COSM PSF Generator" "EstimationGui;COSM Estimation" "ToolsGui;COSM Tools")

INCLUDE(CPack)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(COSM_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
ADD_DEFINITIONS(-Wall -DCOSM_VERSION="${COSM_VERSION}")

SUBDIRS(util psf est tools viewer)
