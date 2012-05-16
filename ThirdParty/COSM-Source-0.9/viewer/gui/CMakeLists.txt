# The name of our project is "VIEWER_GUI". CMakeLists files in this
# project can refer to the root source directory of the project as 
# ${VIEWER_GUI_SOURCE_DIR} and to the root binary directory of the
# project as ${VIEWER_GUI_BINARY_DIR}.
#
PROJECT (VIEWER_GUI)

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

FIND_PATH( BLITZ_DIR $ENV{BLITZ_DIR} [DOC "Blitz directory path"])

INCLUDE_DIRECTORIES (
    ${VIEWER_GUI_SOURCE_DIR} 
    ${VIEWER_GUI_SOURCE_DIR}/.. 
    ${VIEWER_GUI_SOURCE_DIR}/../../util
    ${VIEWER_GUI_SOURCE_DIR}/../../util/fltk
    ${VIEWER_GUI_SOURCE_DIR}/../../util/vtk
    ${FLTK_INCLUDE_DIR}
    ${BLITZ_DIR}
    ${VTK_INCLUDE_DIR}
)

LINK_DIRECTORIES ( 
    ${VIEWER_GUI_BINARY_DIR}/../../lib 
)

SET (VIEWER_GUI_FLUID_SRCS ViewerGUI.fl ViewerAbout.fl)
FLTK_WRAP_UI (CosmViewer ${VIEWER_GUI_FLUID_SRCS})
SET (VIEWER_GUI_SRCS vtkWUReader vtkWUWriter ViewerBase ViewerMain )
IF(UNIX)
    SET(FLTK_EXTRA fltk_jpeg fltk_png fltk_z)
ENDIF(UNIX)

# Add the target
SET(EXECUTABLE_OUTPUT_PATH ${VIEWER_GUI_BINARY_DIR}/../../bin)
ADD_EXECUTABLE (CosmViewer ${CosmViewer_FLTK_UI_SRCS} ${VIEWER_GUI_SRCS})
TARGET_LINK_LIBRARIES (CosmViewer wuheader fltk_contrib vtk_fltk ${FLTK_LIBRARIES} ${FLTK_EXTRA} vtkWidgets vtkHybrid vtkIO )

# install executable
INSTALL_TARGETS(/bin CosmViewer)
