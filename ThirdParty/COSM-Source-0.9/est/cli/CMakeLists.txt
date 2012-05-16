# The name of our project is "EST_CLI". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${EST_CLI_SOURCE_DIR} and to the root binary directory of the
# project as ${EST_CLI_BINARY_DIR}.
#
PROJECT(EST_CLI)
 
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])
FIND_PATH( FFTW_DIR $ENV{FFTW_DIR} [DOC "fftw directory path"])
FIND_PATH( TCLAP_DIR $ENV{TCLAP_DIR} [DOC "TCLAP directory path"])

ADD_DEFINITIONS(-ftemplate-depth-30)

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
    ${EST_CLI_SOURCE_DIR} 
    ${EST_CLI_SOURCE_DIR}/.. 
    ${EST_CLI_SOURCE_DIR}/../../util 
    ${EST_CLI_SOURCE_DIR}/../../util/blitz 
    ${BLITZ_DIR}
    ${TCLAP_DIR}/include
    ${FFTW_INCLUDE_DIR}
)

LINK_DIRECTORIES ( 
    ${EST_CLI_BINARY_DIR}/../../lib
    ${FFTW_LIBRARY_DIR}
)
 
SET(EXECUTABLE_OUTPUT_PATH ${EST_CLI_BINARY_DIR}/../../bin)
 
# source files for cli
SET(EST_CLI_SRCS estMain)
 
# create executable for cli
ADD_EXECUTABLE (EstimationCli ${EST_CLI_SRCS})
TARGET_LINK_LIBRARIES (EstimationCli wuheader ${FFTW_LIBRARY})

# install executable
INSTALL_TARGETS(/bin EstimationCli)
