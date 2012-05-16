# The name of our project is "FLTK_CONTRIB". CMakeLists files in this
# project can refer to the root source directory of the project as
# ${FLTK_CONTRIB_SOURCE_DIR} and to the root binary directory of the
# project as ${FLTK_CONTRIB_BINARY_DIR}.
#
PROJECT(VTK_FLTK)
 
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

FIND_PACKAGE(VTK)
IF(NOT VTK_DIR)
  MESSAGE(FATAL_ERROR "Please set VTK_DIR.")
ENDIF(NOT VTK_DIR)
INCLUDE(${VTK_USE_FILE})

INCLUDE_DIRECTORIES ( 
    ${VTK_FLTK_SOURCE_DIR} 
    ${VTK_FLTK_SOURCE_DIR}/.. 
    ${VTK_FLTK_INCLUDE_DIRS} 
    ${FLTK_INCLUDE_DIR}
    ${VTK_INCLUDE_DIR}
)
 
SET(LIBRARY_OUTPUT_PATH ${VTK_FLTK_BINARY_DIR}/../../lib)
 
# source files for vtk fltk
SET(VTK_FLTK_SRCS 
    Fl_VTK_Window.cxx                       
    vtkFLTKObjectFactory.cxx                
    vtkFLTKOpenGLRenderWindow.cxx
    vtkFLTKRenderWindowInteractor.cxx
)
 
# create library for wuheader
ADD_LIBRARY(vtk_fltk ${VTK_FLTK_SRCS})
