# The name of our project is "FLTK_CONTRIB". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${FLTK_CONTRIB_SOURCE_DIR} and to the root binary directory of the
# project as ${FLTK_CONTRIB_BINARY_DIR}.
#
PROJECT(FLTK_CONTRIB)
 
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

INCLUDE_DIRECTORIES ( ${FLTK_CONTRIB_SOURCE_DIR} ${FLTK_CONTRIB_SOURCE_DIR}/.. ${FLTK_INCLUDE_DIR} )
 
SET(LIBRARY_OUTPUT_PATH ${FLTK_CONTRIB_BINARY_DIR}/../../lib)
 
# source files for fltk contrib

SET(FLTK_CONTRIB_SRCS 
    Table
    Cartesian 

)
 
# create library for wuheader
ADD_LIBRARY(fltk_contrib ${FLTK_CONTRIB_SRCS})
