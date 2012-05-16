project(ITKMicroscopyPSFToolkit)

set( COSM_PATH
  ${ITKMicroscopyPSFToolkit_SOURCE_DIR}/ThirdParty/COSM-Source-0.9
  )

set(ITKMicroscopyPSFToolkit_INCLUDE_DIRS
  "${COSM_PATH}/psf"
  "${COSM_PATH}/util"
)

itk_module_impl()
