# The name of our project is "QUADPACK". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${QUADPACK_SOURCE_DIR} and to the root binary directory of the
# project as ${QUADPACK_BINARY_DIR}.
#
PROJECT(UTILS)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE debug CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SUBDIRS(quadpack tinyxml wu fltk vtk)
