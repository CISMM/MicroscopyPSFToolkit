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

#include "itkHaeberleCOSMOSPointSpreadFunctionImageSource.h"

#include "itkImageFileWriter.h"
#include "itkTestingMacros.h"

#include <cstdlib>

int itkHaeberleCOSMOSPointSpreadFunctionImageSourceTest(int argc, char * argv[])
{
  if ( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " <output file name>" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image< double, 3 >                                        ImageType;
  typedef itk::HaeberleCOSMOSPointSpreadFunctionImageSource< ImageType > SourceType;

  SourceType::Pointer source = SourceType::New();
  SourceType::SizeType size = {{32, 32, 16}};
  source->SetSize( size );

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[1] );
  writer->SetInput( source->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
