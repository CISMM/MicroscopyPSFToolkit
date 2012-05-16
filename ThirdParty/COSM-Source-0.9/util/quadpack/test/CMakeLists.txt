# The name of our project is "QUADPACK_TEST". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${QUADPACK_TEST_SOURCE_DIR} and to the root binary directory of the
# project as ${QUADPACK_TEST_BINARY_DIR}.
#
PROJECT(QUADPACK_TEST)
 
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE debug CACHE STRING
      "Choose the type of build, options are: debug release"
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF (NOT WIN32)
    SET( MATH m )
ENDIF(NOT WIN32)

INCLUDE_DIRECTORIES (
    ${QUADPACK_TEST_SOURCE_DIR} 
    ${QUADPACK_TEST_SOURCE_DIR}/.. 
)

LINK_DIRECTORIES ( 
    ${QUADPACK_TEST_BINARY_DIR}/../../../lib
)
 
SET(OBJECT_OUTPUT_PATH ${CMAKE_BUILD_TYPE})

# create executable for quadpack test
ADD_EXECUTABLE (dq4tst dq4tst)
TARGET_LINK_LIBRARIES (dq4tst quadpack ${MATH} )

ADD_EXECUTABLE (dqaget dqaget)
TARGET_LINK_LIBRARIES (dqaget quadpack ${MATH})

ADD_EXECUTABLE (dqagit dqagit)
TARGET_LINK_LIBRARIES (dqagit quadpack ${MATH})

ADD_EXECUTABLE (dqagpt dqagpt)
TARGET_LINK_LIBRARIES (dqagpt quadpack ${MATH})

ADD_EXECUTABLE (dqagst dqagst)
TARGET_LINK_LIBRARIES (dqagst quadpack ${MATH})

ADD_EXECUTABLE (dqagtst dqagtst)
TARGET_LINK_LIBRARIES (dqagtst quadpack ${MATH})

ADD_EXECUTABLE (dqawct dqawct)
TARGET_LINK_LIBRARIES (dqawct quadpack ${MATH})

ADD_EXECUTABLE (dqawft dqawft)
TARGET_LINK_LIBRARIES (dqawft quadpack ${MATH})

ADD_EXECUTABLE (dqawot dqawot)
TARGET_LINK_LIBRARIES (dqawot quadpack ${MATH})

ADD_EXECUTABLE (dqawst dqawst)
TARGET_LINK_LIBRARIES (dqawst quadpack ${MATH})

ADD_EXECUTABLE (dqct dqct)
TARGET_LINK_LIBRARIES (dqct quadpack ${MATH})

ADD_EXECUTABLE (dqngt dqngt)
TARGET_LINK_LIBRARIES (dqngt quadpack ${MATH})
