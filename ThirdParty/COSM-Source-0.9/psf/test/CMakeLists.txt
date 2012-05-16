# The name of our project is "PSF_TEST". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${PSF_TEST_SOURCE_DIR} and to the root binary directory of the
# project as ${PSF_TEST_BINARY_DIR}.
#
PROJECT(PSF_TEST)
 
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])

ADD_DEFINITIONS(-DTINYXML_USE_STL -Wall)

INCLUDE_DIRECTORIES (
    ${PSF_TEST_SOURCE_DIR} 
    ${PSF_TEST_SOURCE_DIR}/.. 
    ${PSF_TEST_SOURCE_DIR}/../../util 
    ${BLITZ_DIR}
)
LINK_DIRECTORIES ( 
    ${PSF_TEST_BINARY_DIR}/../../lib
)
 
# source files for test
SET(PSF_TEST_SRCS psfMain)
 
# create executable for test
ADD_EXECUTABLE (besselMain besselMain)
TARGET_LINK_LIBRARIES (besselMain psf)

ADD_EXECUTABLE (integratorMain integratorMain)
TARGET_LINK_LIBRARIES (integratorMain psf quadpack)

#ADD_EXECUTABLE (interpolatorMain interpolatorMain)
#TARGET_LINK_LIBRARIES (interpolatorMain psf)

ADD_EXECUTABLE (opdMain opdMain)
TARGET_LINK_LIBRARIES (opdMain psf quadpack)

#ADD_EXECUTABLE (timeMain timeMain)
#TARGET_LINK_LIBRARIES (timeMain psf quadpack)

