/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkGibsonLanniPSFImageSource.cxx,v $
  Language:  C++
  Date:      $Date: 2010/05/17 15:41:35 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
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
