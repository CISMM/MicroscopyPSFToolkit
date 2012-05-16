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
#ifndef __itkCOSMOSPointSpreadFunctionImageSource_hxx
#define __itkCOSMOSPointSpreadFunctionImageSource_hxx

#include "itkCOSMOSPointSpreadFunctionImageSource.h"

namespace itk {

template< class TOutputImage >
COSMOSPointSpreadFunctionImageSource< TOutputImage >
::COSMOSPointSpreadFunctionImageSource()
{
  m_Magnification                     = 60;    // unitless
  m_NumericalAperture                 = 1.4;   // unitless
  m_DesignCoverSlipRefractiveIndex    = 1.522; // unitless
  m_ActualCoverSlipRefractiveIndex    = 1.522; // unitless
  m_DesignCoverSlipThickness          = 170.0; // in micrometers
  m_ActualCoverSlipThickness          = 170.0; // in micrometers
  m_DesignImmersionOilRefractiveIndex = 1.515; // unitless
  m_ActualImmersionOilRefractiveIndex = 1.515; // unitless
  m_DesignImmersionOilThickness       = 100.0; // in micrometers

  m_DesignSpecimenLayerRefractiveIndex    =  1.33; // unitless
  m_ActualSpecimenLayerRefractiveIndex    =  1.33; // unitless
  m_ActualPointSourceDepthInSpecimenLayer =   0.0; // in micrometers

  m_ShearX = 0.0; // nm in X vs. nm in Z
  m_ShearY = 0.0; // nm in Y vs. nm in Z
}

template< class TOutputImage >
COSMOSPointSpreadFunctionImageSource< TOutputImage >
::~COSMOSPointSpreadFunctionImageSource()
{
}

template< class TOutputImage >
void
COSMOSPointSpreadFunctionImageSource< TOutputImage >
::SetParameters(const ParametersType& parameters)
{
  int index = 0;

  this->SetMagnification(parameters[index++]);
  this->SetNumericalAperture(parameters[index++]);
  this->SetEmissionWavelength(parameters[index++]);

  this->SetDesignCoverSlipRefractiveIndex(parameters[index++]);
  this->SetActualCoverSlipRefractiveIndex(parameters[index++]);
  this->SetDesignCoverSlipThickness(parameters[index++]);
  this->SetActualCoverSlipThickness(parameters[index++]);
  this->SetDesignImmersionOilRefractiveIndex(parameters[index++]);
  this->SetActualImmersionOilRefractiveIndex(parameters[index++]);
  this->SetDesignImmersionOilThickness(parameters[index++]);

  this->SetDesignSpecimenLayerRefractiveIndex(parameters[index++]);
  this->SetActualSpecimenLayerRefractiveIndex(parameters[index++]);
  this->SetActualPointSourceDepthInSpecimenLayer(parameters[index++]);

  this->SetShearX(parameters[index++]);
  this->SetShearY(parameters[index++]);
}

template< class TOutputImage >
typename COSMOSPointSpreadFunctionImageSource< TOutputImage >::ParametersType
COSMOSPointSpreadFunctionImageSource< TOutputImage >
::GetParameters() const
{
  ParametersType parameters( this->GetNumberOfParameters() );

  int index = 0;

  parameters[index++] = this->GetMagnification();
  parameters[index++] = this->GetNumericalAperture();
  parameters[index++] = this->GetEmissionWavelength();

  parameters[index++] = this->GetDesignCoverSlipRefractiveIndex();
  parameters[index++] = this->GetActualCoverSlipRefractiveIndex();
  parameters[index++] = this->GetDesignCoverSlipThickness();
  parameters[index++] = this->GetActualCoverSlipThickness();
  parameters[index++] = this->GetDesignImmersionOilRefractiveIndex();
  parameters[index++] = this->GetActualImmersionOilRefractiveIndex();
  parameters[index++] = this->GetDesignImmersionOilThickness();

  parameters[index++] = this->GetDesignSpecimenLayerRefractiveIndex();
  parameters[index++] = this->GetActualSpecimenLayerRefractiveIndex();
  parameters[index++] = this->GetActualPointSourceDepthInSpecimenLayer();

  parameters[index++] = this->GetShearX();
  parameters[index++] = this->GetShearY();

  return parameters;
}

template< class TOutputImage >
unsigned int
COSMOSPointSpreadFunctionImageSource< TOutputImage >
::GetNumberOfParameters() const
{
  return 15;
}

template< class TOutputImage >
void
COSMOSPointSpreadFunctionImageSource< TOutputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  os << indent << "Magnification: " << m_Magnification << "\n";
  os << indent << "NumericalAperture: " << m_NumericalAperture << "\n";
  os << indent << "EmissionWavelength: " << m_EmissionWavelength << "\n";
  os << indent << "DesignCoverSlipRefractiveIndex: " << m_DesignCoverSlipRefractiveIndex << "\n";
  os << indent << "ActualCoverSlipRefractiveIndex: " << m_ActualSpecimenLayerRefractiveIndex << "\n";
  os << indent << "DesignCoverSlipThickness: " << m_DesignCoverSlipThickness << "\n";
  os << indent << "ActualCoverSlipThickness: " << m_ActualCoverSlipThickness << "\n";
  os << indent << "DesignImmersionOilRefractiveIndex: " << m_DesignImmersionOilRefractiveIndex << "\n";
  os << indent << "ActualImmersionOilRefractiveIndex: " << m_ActualImmersionOilRefractiveIndex << "\n";
  os << indent << "DesignImmersionOilThickness: " << m_DesignImmersionOilThickness << "\n";
  os << indent << "DesignSpecimenLayerRefractiveIndex: " << m_DesignSpecimenLayerRefractiveIndex << "\n";
  os << indent << "ActualSpecimenLayerRefractiveIndex: " << m_ActualSpecimenLayerRefractiveIndex << "\n";
  os << indent << "ActualPointSourceDepthInSpecimenLayer: " << m_ActualPointSourceDepthInSpecimenLayer << "\n";
  os << indent << "ShearX: " << m_ShearX << "\n";
  os << indent << "ShearY: " << m_ShearY << "\n";
}


} // end namespace itk

#endif

