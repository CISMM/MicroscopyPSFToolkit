# The name of our project is "PSF_CLI". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${PSF_CLI_SOURCE_DIR} and to the root binary directory of the
# project as ${PSF_CLI_BINARY_DIR}.
#
PROJECT(PSF_CLI)
 
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# Include the FLTK CMake configuration.
IF (NOT FIND_FLTK_FILE)
  FIND_FILE (FIND_FLTK_FILE FindFLTK.cmake
    ${CMAKE_ROOT}/Modules
    )
ENDIF (NOT FIND_FLTK_FILE)

IF (FIND_FLTK_FILE)
  INCLUDE (${FIND_FLTK_FILE})
ELSE (FIND_FLTK_FILE)
  MESSAGE(SEND_ERROR "FindFLTK.cmake not found!")
ENDIF (FIND_FLTK_FILE)

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz++ directory path"])
FIND_PATH( TCLAP_DIR $ENV{TCLAP_DIR} [DOC "TCLAP directory path"])


ADD_DEFINITIONS(-DTINYXML_USE_STL)

INCLUDE_DIRECTORIES (
    ${PSF_CLI_SOURCE_DIR} 
    ${PSF_CLI_SOURCE_DIR}/.. 
    ${PSF_CLI_SOURCE_DIR}/../../util 
    ${FLTK_INCLUDE_DIR}
    ${BLITZ_DIR}
    ${TCLAP_DIR}/include
)

LINK_DIRECTORIES ( 
    ${PSF_CLI_BINARY_DIR}/../../lib
)
 
SET(EXECUTABLE_OUTPUT_PATH ${PSF_CLI_BINARY_DIR}/../../bin)
 
# source files for cli
SET(PSF_CLI_SRCS psfMain)
 
# create executable for cli
ADD_EXECUTABLE (PsfGeneratorCli ${PSF_CLI_SRCS})
TARGET_LINK_LIBRARIES (PsfGeneratorCli psf wuheader quadpack tinyxml)

# install executable
INSTALL_TARGETS(/bin PsfGeneratorCli)

