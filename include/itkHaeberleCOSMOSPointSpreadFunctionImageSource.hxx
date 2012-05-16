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

#ifndef __itkHaeberleCOSMOSPointSpreadFunctionImageSource_hxx
#define __itkHaeberleCOSMOSPointSpreadFunctionImageSource_hxx

#include "itkHaeberleCOSMOSPointSpreadFunctionImageSource.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkMath.h"
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"

#include <algorithm>


namespace itk
{

//----------------------------------------------------------------------------
template <class TOutputImage>
HaeberleCOSMOSPointSpreadFunctionImageSource<TOutputImage>
::HaeberleCOSMOSPointSpreadFunctionImageSource()
{
}


//----------------------------------------------------------------------------
template <class TOutputImage>
HaeberleCOSMOSPointSpreadFunctionImageSource<TOutputImage>
::~HaeberleCOSMOSPointSpreadFunctionImageSource()
{
}

//----------------------------------------------------------------------------
template < class TOutputImage >
void
HaeberleCOSMOSPointSpreadFunctionImageSource<TOutputImage>
::GenerateData()
{
  m_HaeberleFunctor = new cosm::HaeberlePsfFunctor< double >(
    1e-3*this->GetActualPointSourceDepthInSpecimenLayer(),
    1e-3*this->GetDesignImmersionOilThickness(),
    1e-3*this->GetDesignImmersionOilThickness(), // I didn't think this
                                              // needed to be defined...
    1e-3*this->GetDesignCoverSlipThickness(),
    1e-3*this->GetActualCoverSlipThickness(),
    this->GetActualSpecimenLayerRefractiveIndex(),
    this->GetDesignImmersionOilRefractiveIndex(),
    this->GetActualImmersionOilRefractiveIndex(),
    this->GetDesignCoverSlipRefractiveIndex(),
    this->GetActualCoverSlipRefractiveIndex(),
    160.0, /*design tube length*/
    160.0, /*actual tube length*/
    this->GetMagnification(),
    this->GetNumericalAperture(),
    1e-6*this->GetEmissionWavelength(),
    1e-6);

  // allocate the output buffer
  OutputImageType * output = this->GetOutput();
  output->SetBufferedRegion( output->GetRequestedRegion() );
  output->Allocate();

  OutputImageRegionType outputRegion = output->GetRequestedRegion();
  ImageRegionIteratorWithIndex< OutputImageType > it(output, outputRegion);

  while ( !it.IsAtEnd() )
    {
    OutputImageIndexType index = it.GetIndex();
    OutputImagePointType point;
    output->TransformIndexToPhysicalPoint( index, point );

    it.Set( this->ComputeSampleValue( point ) );
    ++it;
    }

  delete m_HaeberleFunctor;
}

//----------------------------------------------------------------------------
template <class TOutputImage>
double
HaeberleCOSMOSPointSpreadFunctionImageSource<TOutputImage>
::ComputeSampleValue(const OutputImagePointType& point)
{
  // Convert from nanometers to millimeters
  OutputImagePointValueType px = point[0] * 1e-6;
  OutputImagePointValueType py = point[1] * 1e-6;
  OutputImagePointValueType pz = point[2] * 1e-6;

  double r = sqrt( (px*px) + (py*py) );
  return static_cast< OutputImagePixelType >( norm( (*m_HaeberleFunctor)( pz, r ) ) );
}


} // end namespace itk

#endif
