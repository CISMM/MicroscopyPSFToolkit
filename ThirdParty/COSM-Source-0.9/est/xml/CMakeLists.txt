# The name of our project is "EST_XML". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${EST_XML_SOURCE_DIR} and to the root binary directory of the
# project as ${EST_XML_BINARY_DIR}.
#
PROJECT(EST_XML)
 
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH( TCLAP_DIR $ENV{TCLAP_DIR} [DOC "TCLAP directory path"])
FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])

ADD_DEFINITIONS(-ftemplate-depth-30)

INCLUDE_DIRECTORIES (
    ${EST_XML_SOURCE_DIR} 
    ${EST_XML_SOURCE_DIR}/.. 
    ${EST_XML_SOURCE_DIR}/../../util 
    ${BLITZ_DIR}
    ${TCLAP_DIR}/include
)

LINK_DIRECTORIES ( 
    ${EST_XML_BINARY_DIR}/../../lib
)
 
SET(EXECUTABLE_OUTPUT_PATH ${EST_XML_BINARY_DIR}/../../bin)
 
# source files for xml
SET(EST_XML_SRCS estXmlMain)
 
# create executable for xml
ADD_EXECUTABLE (EstXmlMain ${EST_XML_SRCS})
TARGET_LINK_LIBRARIES (EstXmlMain tinyxml)
