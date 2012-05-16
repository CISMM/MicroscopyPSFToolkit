# The name of our project is "TOOLS_CLI". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${TOOLS_CLI_SOURCE_DIR} and to the root binary directory of the
# project as ${TOOLS_CLI_BINARY_DIR}.
#
PROJECT(TOOLS_CLI)
 
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

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])
FIND_PATH( TCLAP_DIR $ENV{TCLAP_DIR} [DOC "TCLAP directory path"])
FIND_PATH( FFTW_DIR $ENV{FFTW_DIR} [DOC "FFTW directory path"])

ADD_DEFINITIONS(-DTINYXML_USE_STL)

IF (NOT WIN32)
    SET( FFTW_LIBRARY fftw3f fftw3 fftw3l)
    SET( FFTW_LIBRARY_DIR ${FFTW_DIR}/lib )
    SET( FFTW_INCLUDE_DIR ${FFTW_DIR}/include )
ELSE (NOT WIN32)
    SET( FFTW_LIBRARY libfftw3f-3 libfftw3-3 libfftw3l-3)
    SET( FFTW_LIBRARY_DIR ${FFTW_DIR} )
    SET( FFTW_INCLUDE_DIR ${FFTW_DIR} )
ENDIF (NOT WIN32)

INCLUDE_DIRECTORIES (
    ${TOOLS_CLI_SOURCE_DIR} 
    ${TOOLS_CLI_SOURCE_DIR}/.. 
    ${TOOLS_CLI_SOURCE_DIR}/../src 
    ${TOOLS_CLI_SOURCE_DIR}/../../est 
    ${TOOLS_CLI_SOURCE_DIR}/../../util 
    ${FLTK_INCLUDE_DIR}
    ${BLITZ_DIR}
    ${TCLAP_DIR}/include
    ${FFTW_INCLUDE_DIR}
)
LINK_DIRECTORIES ( 
    ${TOOLS_CLI_BINARY_DIR}/../../lib
    ${FFTW_LIBRARY_DIR}
)
 
SET(EXECUTABLE_OUTPUT_PATH ${TOOLS_CLI_BINARY_DIR}/../../bin)
 
# create executable for tool clii

ADD_EXECUTABLE (ToolsCli toolsMain)
TARGET_LINK_LIBRARIES (ToolsCli tools wuheader ${FFTW_LIBRARY})

# install executable
INSTALL_TARGETS(/bin ToolsCli)
