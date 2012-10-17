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

#ifndef __itkGibsonLanniCOSMOSPointSpreadFunctionImageSource_hxx
#define __itkGibsonLanniCOSMOSPointSpreadFunctionImageSource_hxx

#include "itkGibsonLanniCOSMOSPointSpreadFunctionImageSource.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkMath.h"
#include "itkObjectFactory.h"
#include "itkProgressReporter.h"

namespace itk
{

//----------------------------------------------------------------------------
template< class TOutputImage >
GibsonLanniCOSMOSPointSpreadFunctionImageSource<TOutputImage>
::GibsonLanniCOSMOSPointSpreadFunctionImageSource()
{
}


//----------------------------------------------------------------------------
template< class TOutputImage >
GibsonLanniCOSMOSPointSpreadFunctionImageSource<TOutputImage>
::~GibsonLanniCOSMOSPointSpreadFunctionImageSource()
{
}

//----------------------------------------------------------------------------
template< class TOutputImage >
void
GibsonLanniCOSMOSPointSpreadFunctionImageSource<TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                       ThreadIdType itkNotUsed(threadId))
{
  cosm::GibsonLaniPsfFunctor< double > gibsonLanniFunctor(
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

  OutputImageType * output = this->GetOutput();

  ImageRegionIteratorWithIndex< OutputImageType > it(output, outputRegionForThread);

  while ( !it.IsAtEnd() )
    {
    OutputImageIndexType index = it.GetIndex();
    OutputImagePointType point;
    output->TransformIndexToPhysicalPoint( index, point );

    // Convert from nanometers to millimeters
    OutputImagePixelType pz = point[2] * 1e-6;
    OutputImagePixelType px = point[0] * 1e-6 + (pz * this->GetShearX());
    OutputImagePixelType py = point[1] * 1e-6 + (pz * this->GetShearY());

    double r = sqrt( (px*px) + (py*py) );
    it.Set( static_cast< OutputImagePixelType >( norm( gibsonLanniFunctor( pz, r ) ) ) );

    ++it;
    }
}

} // end namespace itk

#endif
