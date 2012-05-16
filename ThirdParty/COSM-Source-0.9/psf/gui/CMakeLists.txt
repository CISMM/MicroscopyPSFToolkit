# The name of our project is "PSF_GUI". CMakeLists files in this
# project can refer to the root source directory of the project as 
# ${PSF_GUI_SOURCE_DIR} and to the root binary directory of the
# project as ${PSF_GUI_BINARY_DIR}.
#
PROJECT (PSF_GUI)

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

ADD_DEFINITIONS(-ftemplate-depth-30 -DNOMINMAX -D_USE_MATH_DEFINES)

INCLUDE_DIRECTORIES (
    ${PSF_GUI_SOURCE_DIR}
    ${PSF_GUI_SOURCE_DIR}/..
    ${PSF_GUI_SOURCE_DIR}/../../util  
    ${PSF_GUI_SOURCE_DIR}/../../util/fltk  
    ${FLTK_INCLUDE_DIR}
    ${BLITZ_DIR}
    ${TCLAP_DIR}/include
)

LINK_DIRECTORIES ( 
    ${PSF_GUI_BINARY_DIR}/../../lib 
)

SET (PSF_GUI_FLUID_SRCS psfGUI.fl psfAbout.fl)
SET (PSF_GUI_SRCS psfMain)
IF(WIN32)
SET (FLTK_EXTRA fltkjpeg fltkpng zlib)
ELSE(WIN32)
SET (FLTK_EXTRA fltk_jpeg fltk_png fltk_z)
ENDIF(WIN32)

# Add the target
FLTK_WRAP_UI (CosmPsf ${PSF_GUI_FLUID_SRCS})
SET(EXECUTABLE_OUTPUT_PATH ${PSF_GUI_BINARY_DIR}/../../bin)
ADD_EXECUTABLE (CosmPsf ${CosmPsf_FLTK_UI_SRCS} ${PSF_GUI_SRCS})
TARGET_LINK_LIBRARIES (CosmPsf psf xcosm wuheader fltk_contrib quadpack tinyxml ${FLTK_LIBRARIES} ${FLTK_EXTRA})

# install executable
INSTALL_TARGETS(/bin CosmPsf)
