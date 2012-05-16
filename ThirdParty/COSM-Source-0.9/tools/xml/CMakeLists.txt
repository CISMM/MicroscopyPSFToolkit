# The name of our project is "TOOLS_XML". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${TOOLS_XML_SOURCE_DIR} and to the root binary directory of the
# project as ${TOOLS_XML_BINARY_DIR}.
#
PROJECT(TOOLS_XML)
 
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH( TCLAP_DIR $ENV{TCLAP_DIR} [DOC "TCLAP directory path"])
FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])

INCLUDE_DIRECTORIES (
    ${TOOLS_XML_SOURCE_DIR} 
    ${TOOLS_XML_SOURCE_DIR}/.. 
    ${TOOLS_XML_SOURCE_DIR}/../../util 
    ${BLITZ_DIR}
    ${TCLAP_DIR}/include
)

LINK_DIRECTORIES ( 
    ${TOOLS_XML_BINARY_DIR}/../../lib
)
 
SET(EXECUTABLE_OUTPUT_PATH ${TOOLS_XML_BINARY_DIR}/../../bin)
 
# source files for xml
SET(TOOLS_XML_SRCS toolsXmlMain)
 
# create executable for xml
ADD_EXECUTABLE (ToolsXmlMain ${TOOLS_XML_SRCS})
TARGET_LINK_LIBRARIES (ToolsXmlMain tinyxml)
