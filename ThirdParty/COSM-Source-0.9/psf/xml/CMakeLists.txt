# The name of our project is "PSF_XML". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${PSF_XML_SOURCE_DIR} and to the root binary directory of the
# project as ${PSF_XML_BINARY_DIR}.
#
PROJECT(PSF_XML)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "BLITZ directory path"])
FIND_PATH( TCLAP_DIR $ENV{TCLAP_DIR} [DOC "TCLAP directory path"])

ADD_DEFINITIONS(-DTIXML_USE_STL)

INCLUDE_DIRECTORIES (
    ${PSF_XML_SOURCE_DIR} 
    ${PSF_XML_SOURCE_DIR}/.. 
    ${PSF_XML_SOURCE_DIR}/../../util 
    ${BLITZ_DIR} 
    ${TCLAP_DIR}/include
)
LINK_DIRECTORIES ( 
    ${PSF_XML_BINARY_DIR}/../../lib 
)
 
# create executable for psf xml test
ADD_EXECUTABLE (XmlMain xmlMain)
TARGET_LINK_LIBRARIES (XmlMain tinyxml)
