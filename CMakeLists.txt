project(ITKMicroscopyPSFToolkit)

message( WARNING "ATTENTION: You have enabled the use of the "
"MicroscopyPSFToolkit. This module uses the COSMOS software which is "
"distributed under a GPL license. By enabling this option, the binary of the ITK "
"libraries that you are going to build will be covered by a GPL license, and so "
"it will be any executable that you link against these libraries." )

set( COSM_PATH
  ${ITKMicroscopyPSFToolkit_SOURCE_DIR}/ThirdParty/COSM-Source-0.9
  )

set(ITKMicroscopyPSFToolkit_INCLUDE_DIRS
  "${COSM_PATH}/psf"
  "${COSM_PATH}/util"
)

itk_module_impl()
