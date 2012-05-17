set( QUADPACK_PATH
  ${COSM_PATH}/util/quadpack
  )
file(GLOB QUADPACK_C_SRC
  ${QUADPACK_PATH}/*.c
)
file(GLOB QUADPACK_CXX_SRC
  ${QUADPACK_PATH}/*.cxx
)
set( QUADPACK_SRC
  ${QUADPACK_C_SRC}
  ${QUADPACK_CXX_SRC}
)

add_library(ITKMicroscopyPSFToolkit ${QUADPACK_SRC})
itk_module_target(ITKMicroscopyPSFToolkit)
