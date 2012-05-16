# The name of our project is "TOOLS_LIB". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${TOOLS_LIB_SOURCE_DIR} and to the root binary directory of the
# project as ${TOOLS_LIB_BINARY_DIR}.
#
PROJECT(TOOLS_LIB)
 
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
FIND_PATH( FFTW_DIR $ENV{FFTW_DIR} [DOC "FFTW directory path"])

IF (NOT WIN32)
    SET( FFTW_LIBRARY fftw3f fftw3 fftw3l)
    SET( FFTW_LIBRARY_DIR ${FFTW_DIR}/.libs )
    SET( FFTW_INCLUDE_DIR ${FFTW_DIR}/api )
ELSE (NOT WIN32)
    SET( FFTW_LIBRARY libfftw3f-3 libfftw3-3 libfftw3l-3)
    SET( FFTW_LIBRARY_DIR ${FFTW_DIR} )
    SET( FFTW_INCLUDE_DIR ${FFTW_DIR} )
ENDIF (NOT WIN32)

ADD_DEFINITIONS(-DTINYXML_USE_STL)

INCLUDE_DIRECTORIES (
    ${TOOLS_LIB_SOURCE_DIR} 
    ${TOOLS_LIB_SOURCE_DIR}/.. 
    ${TOOLS_LIB_SOURCE_DIR}/../../est 
    ${TOOLS_LIB_SOURCE_DIR}/../../util 
    ${BLITZ_DIR}
    ${FFTW_INCLUDE_DIR}
)

SET (TOOLS_LIB_SRCS toolMgr )

# create library for tools
SET(LIBRARY_OUTPUT_PATH ${TOOLS_LIB_BINARY_DIR}/../../lib)
ADD_LIBRARY(tools ${TOOLS_LIB_SRCS})
