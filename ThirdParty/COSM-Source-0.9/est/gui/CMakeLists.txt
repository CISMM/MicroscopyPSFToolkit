# The name of our project is "EST_GUI". CMakeLists files in this
# project can refer to the root source directory of the project as 
# ${EST_GUI_SOURCE_DIR} and to the root binary directory of the
# project as ${EST_GUI_BINARY_DIR}.
#
PROJECT (EST_GUI)

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

ADD_DEFINITIONS(-ftemplate-depth-30 -DNOMINMAX )

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
    ${EST_GUI_SOURCE_DIR} 
    ${EST_GUI_SOURCE_DIR}/..
    ${EST_GUI_SOURCE_DIR}/../../util  
    ${EST_GUI_SOURCE_DIR}/../../util/blitz  
    ${EST_GUI_SOURCE_DIR}/../../util/fltk  
    ${FFTW_INCLUDE_DIR}
    ${FLTK_INCLUDE_DIR}
    ${BLITZ_DIR}
    ${TCLAP_DIR}/include
)

LINK_DIRECTORIES ( 
    ${EST_GUI_BINARY_DIR}/../../lib
    ${FFTW_LIBRARY_DIR}
)

SET (EST_GUI_FLUID_SRCS estGUI.fl estAbout.fl)
SET (EST_GUI_SRCS estMain estAlgo)
IF(UNIX)
    SET(FLTK_EXTRA fltk_jpeg fltk_png fltk_z)
ENDIF(UNIX)

# Add the target
FLTK_WRAP_UI (CosmEstimation ${EST_GUI_FLUID_SRCS})
SET(EXECUTABLE_OUTPUT_PATH ${EST_GUI_BINARY_DIR}/../../bin)
ADD_EXECUTABLE (CosmEstimation ${CosmEstimation_FLTK_UI_SRCS} ${EST_GUI_SRCS})
TARGET_LINK_LIBRARIES (CosmEstimation wuheader tinyxml fltk_contrib ${FFTW_LIBRARY} ${FLTK_LIBRARY} ${FLTK_EXTRA})

# install executable
INSTALL_TARGETS(/bin CosmEstimation)
IF(WIN32)
    INSTALL_FILES(/bin FILES ${FFTW_LIBRARY_DIR}/libfftw3f-3.dll ${FFTW_LIBRARY_DIR}/libfftw3-3.dll ${FFTW_LIBRARY_DIR}/libfftw3l-3.dll)
ENDIF(WIN32)
