/****************************************************************************
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 ****************************************************************************/

#include "itkGibsonLanniCOSMOSPointSpreadFunctionImageSource.h"

#include "itkImageFileWriter.h"
#include "itkTestingMacros.h"

#include <cstdlib>

int itkGibsonLanniCOSMOSPointSpreadFunctionImageSourceTest(int argc, char * argv[])
{
  if ( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " <output file name>" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image< double, 3 >                                           ImageType;
  typedef itk::GibsonLanniCOSMOSPointSpreadFunctionImageSource< ImageType > SourceType;

  SourceType::Pointer source = SourceType::New();
  SourceType::SizeType size = {{32, 32, 16}};
  source->SetSize( size );

  SourceType::SpacingType spacing;
  spacing[0] = 65.0;
  spacing[1] = 65.0;
  spacing[2] = 200.0;
  source->SetSpacing( spacing );

  SourceType::PointType origin;
  for (int i = 0; i < ImageType::ImageDimension; ++i)
    {
    origin[i] = -0.5 * ( spacing[i] * static_cast< double >(size[i]-1) );
    }
  source->SetOrigin( origin );

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( source->GetOutput() );
  writer->Update();

  // Now exercise the setters/getters
  // Check that the number of parameters is what we expect
  TEST_SET_GET_VALUE( 15, source->GetNumberOfParameters() );

  // Test setters and getters
  SourceType::ParametersValueType param;
  param = 0.0;
  source->SetMagnification( param );
  TEST_SET_GET_VALUE( param, source->GetMagnification() );

  param = 1.0;
  source->SetNumericalAperture( param );
  TEST_SET_GET_VALUE( param, source->GetNumericalAperture() );

  param = 2.0;
  source->SetEmissionWavelength( param );
  TEST_SET_GET_VALUE( param, source->GetEmissionWavelength() );

  param = 3.0;
  source->SetDesignCoverSlipRefractiveIndex( param );
  TEST_SET_GET_VALUE( param, source->GetDesignCoverSlipRefractiveIndex() );

  param = 4.0;
  source->SetActualCoverSlipRefractiveIndex( param );
  TEST_SET_GET_VALUE( param, source->GetActualCoverSlipRefractiveIndex() );

  param = 5.0;
  source->SetDesignCoverSlipThickness( param );
  TEST_SET_GET_VALUE( param, source->GetDesignCoverSlipThickness() );

  param = 6.0;
  source->SetActualCoverSlipThickness( param );
  TEST_SET_GET_VALUE( param, source->GetActualCoverSlipThickness() );

  param = 7.0;
  source->SetDesignImmersionOilRefractiveIndex( param );
  TEST_SET_GET_VALUE( param, source->GetDesignImmersionOilRefractiveIndex() );

  param = 8.0;
  source->SetActualImmersionOilRefractiveIndex( param );
  TEST_SET_GET_VALUE( param, source->GetActualImmersionOilRefractiveIndex() );

  param = 9.0;
  source->SetDesignImmersionOilThickness( param );
  TEST_SET_GET_VALUE( param, source->GetDesignImmersionOilThickness() );

  param = 10.0;
  source->SetDesignSpecimenLayerRefractiveIndex( param );
  TEST_SET_GET_VALUE( param, source->GetDesignSpecimenLayerRefractiveIndex() );

  param = 11.0;
  source->SetActualSpecimenLayerRefractiveIndex( param );
  TEST_SET_GET_VALUE( param, source->GetActualSpecimenLayerRefractiveIndex() );

  param = 12.0;
  source->SetActualPointSourceDepthInSpecimenLayer( param );
  TEST_SET_GET_VALUE( param, source->GetActualPointSourceDepthInSpecimenLayer() );

  param = 13.0;
  source->SetShearX( param );
  TEST_SET_GET_VALUE( param, source->GetShearX() );

  param = 14.0;
  source->SetShearY( param );
  TEST_SET_GET_VALUE( param, source->GetShearY() );

  // Now check that the parameters retrieved through GetParameters()
  // is what we expect
  SourceType::ParametersType params = source->GetParameters();
  for (unsigned int i = 0; i < source->GetNumberOfParameters(); ++i)
    {
    SourceType::ParametersValueType expectedValue =
      static_cast< SourceType::ParametersValueType >( i );
    if ( params[i] != expectedValue )
      {
      std::cerr << "Expected parameter " << i << " to have value "
                << expectedValue << ", was " << params[i] << std::endl;
      return EXIT_FAILURE;
      }
    }

  // Set these parameters and check that we get them back
  source->SetParameters( params );
  SourceType::ParametersType newParams = source->GetParameters();
  if ( params != newParams )
    {
    std::cerr << "Parameters did not match.\n";
    std::cerr << "Expected " << params << "\n";
    std::cerr << "but got " << newParams << "\n";
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
