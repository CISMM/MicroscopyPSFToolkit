# The name of our project is "WU". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${WU_SOURCE_DIR} and to the root binary directory of the
# project as ${WU_BINARY_DIR}.
#
PROJECT(WU)
 
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])

INCLUDE_DIRECTORIES (${WU_SOURCE_DIR} ${WU_SOURCE_DIR}/.. ${BLITZ_DIR})
 
SET(LIBRARY_OUTPUT_PATH ${WU_BINARY_DIR}/../../lib)
 
# source files for wu
SET(WU_SRCS wuHeader wuImage)
 
# create library for wuheader
ADD_LIBRARY(wuheader ${WU_SRCS})
