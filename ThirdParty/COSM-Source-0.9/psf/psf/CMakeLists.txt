# The name of our project is "PSF_LIB". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${PSF_LIB_SOURCE_DIR} and to the root binary directory of the
# project as ${PSF_LIB_BINARY_DIR}.
#
PROJECT(PSF_LIB)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

INCLUDE_DIRECTORIES (
    ${PSF_LIB_SOURCE_DIR} 
    ${PSF_LIB_SOURCE_DIR}/../../util 
)
 
# source files for psf
SET(PSF_LIB_SRCS 
    dqagFunction
)

# create library for psf
SET(LIBRARY_OUTPUT_PATH ${PSF_LIB_BINARY_DIR}/../../lib)
ADD_LIBRARY(psf ${PSF_LIB_SRCS})
