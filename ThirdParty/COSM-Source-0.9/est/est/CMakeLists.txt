# The name of our project is "EST_LIB". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${EST_LIB_SOURCE_DIR} and to the root binary directory of the
# project as ${EST_LIB_BINARY_DIR}.
#
PROJECT(EST_LIB)
 
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])
FIND_PATH( FFTW_DIR $ENV{FFTW_DIR} [DOC "fftw directory path"])

ADD_DEFINITIONS(-ftemplate-depth-30)

IF (NOT WIN32)
    SET( FFTW_LIBRARY fftw3f fftw3 fftw3l)
    SET( FFTW_LIBRARY_DIR ${FFTW_DIR}/.libs )
    SET( FFTW_INCLUDE_DIR ${FFTW_DIR}/api )
ELSE (NOT WIN32)
    SET( FFTW_LIBRARY libfftw3f-3 libfftw3-3 libfftw3l-3)
    SET( FFTW_LIBRARY_DIR ${FFTW_DIR} )
    SET( FFTW_INCLUDE_DIR ${FFTW_DIR} )
ENDIF (NOT WIN32)

INCLUDE_DIRECTORIES (
    $(BLITZ_DIR)
    ${EST_LIB_SOURCE_DIR} 
    ${EST_LIB_SOURCE_DIR}/.. 
    ${EST_LIB_SOURCE_DIR}/../../util 
    ${EST_LIB_SOURCE_DIR}/../../util/blitz 
    ${BLITZ_DIR}
    ${FFTW_INCLUDE_DIR}
)

LINK_DIRECTORIES ( 
    ${EST_LIB_BINARY_DIR}/../../lib
    ${FFTW_LIBRARY_DIR}
)
 
# source files for est
SET(EST_LIB_SRCS testEstMain)
 
# create executable for est
ADD_EXECUTABLE (TestEstMain ${EST_LIB_SRCS})
TARGET_LINK_LIBRARIES (TestEstMain wuheader ${FFTW_LIBRARY} )
