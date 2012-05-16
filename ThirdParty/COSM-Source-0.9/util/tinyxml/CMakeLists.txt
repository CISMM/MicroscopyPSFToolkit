# The name of our project is "TINYXML". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${TINYXML_SOURCE_DIR} and to the root binary directory of the
# project as ${TINYXML_BINARY_DIR}.
#
PROJECT(TINYXML)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

ADD_DEFINITIONS(-DTIXML_USE_STL)
INCLUDE_DIRECTORIES (
    ${TINYXML_SOURCE_DIR} ${TINYXML_SOURCE_DIR}/.. 
)

LINK_DIRECTORIES ( ${TINYXML_BINARY_DIR}/../../lib )
 
# source files for tinyxml
SET(TINYXML_SRCS 
    tinyxml tinyxmlparser tinyxmlerror tinystr
)

# create library for tinyxml
SET(LIBRARY_OUTPUT_PATH ${TINYXML_BINARY_DIR}/../../lib )
ADD_LIBRARY(tinyxml ${TINYXML_SRCS})
 
# create executable for tinyxml test
ADD_EXECUTABLE (xmltest xmltest)
TARGET_LINK_LIBRARIES (xmltest tinyxml)
