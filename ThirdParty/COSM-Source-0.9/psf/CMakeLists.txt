# The name of our project is "PSF". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${PSF_SOURCE_DIR} and to the root binary directory of the
# project as ${PSF_BINARY_DIR}.
#
PROJECT(PSF)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SUBDIRS(psf xml xcosm gui test)
