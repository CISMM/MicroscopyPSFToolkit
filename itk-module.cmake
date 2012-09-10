set(DOCUMENTATION "This module contains point-spread function models
for microscopy. The models are generated with portions of the COSMOS
software available from http://cirl.memphis.edu/cosmos.php.")

itk_module(ITKMicroscopyPSFToolkit
 DEPENDS
  ITKCommon
  ITKImageSources
 TEST_DEPENDS
  ITKTestKernel
 DESCRIPTION
  "${DOCUMENTATION}"
)
