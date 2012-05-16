# The name of our project is "QUADPACK". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${QUADPACK_SOURCE_DIR} and to the root binary directory of the
# project as ${QUADPACK_BINARY_DIR}.
#
PROJECT(QUADPACK)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE release CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SUBDIRS(test)

INCLUDE_DIRECTORIES (${QUADPACK_SOURCE_DIR} ${QUADPACK_SOURCE_DIR}/.. )
 
# source files for quadpack
SET(QUADPACK_SRCS 
   dqag
   dqage
   dqagi
   dqagp
   dqags
   dqawc
   dqawce
   dqawf
   dqawfe
   dqawo
   dqaws
   dqawse
   dqc25c
   dqc25o
   dqc25s
   dqcheb
   dqext
   dqfour
   dqk15
   dqk15i
   dqk15w
   dqk21
   dqk31
   dqk41
   dqk51
   dqk61
   dqmomo
   dqng
   dqsort
   dqwgt
)

# create library for quadpack
SET(LIBRARY_OUTPUT_PATH ${QUADPACK_BINARY_DIR}/../../lib)
ADD_LIBRARY(quadpack ${QUADPACK_SRCS})
