# The name of our project is "PSF_XCOSM_LIB". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${PSF_XCOSM_LIB_SOURCE_DIR} and to the root binary directory of the
# project as ${PSF_XCOSM_LIB_BINARY_DIR}.
#
PROJECT(PSF_XCOSM_LIB)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])

INCLUDE_DIRECTORIES (
    ${PSF_XCOSM_LIB_SOURCE_DIR} 
    ${PSF_XCOSM_LIB_SOURCE_DIR}/.. 
    ${PSF_XCOSM_LIB_SOURCE_DIR}/../../util 
    ${BLITZ_DIR}
)
 
# source files for psf
SET(PSF_XCOSM_LIB_SRCS 
    rdlinesum
    rdcircsum
    rotnone
    rotsum
    rotdiskline
    rotdicline
    rotinfdicline
    rotdiskcirc
    spline
    real_fft3d
    write_dataset
    misc
    tmp
    fftn
)

# create library for psf
SET(LIBRARY_OUTPUT_PATH ${PSF_XCOSM_LIB_BINARY_DIR}/../../lib)
ADD_LIBRARY(xcosm ${PSF_XCOSM_LIB_SRCS})

LINK_DIRECTORIES ( 
    ${PSF_XCOSM_LIB_BINARY_DIR}/../../lib
)
 
# source files for est
SET(PSF_XCOSM_LIB_SRCS TestXcosmMain)
 
# create executable for est
ADD_EXECUTABLE (TestXcosmMain ${PSF_XCOSM_LIB_SRCS})
TARGET_LINK_LIBRARIES (TestXcosmMain xcosm )
