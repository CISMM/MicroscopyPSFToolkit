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
#ifndef __itkCOSMOSPointSpreadFunctionImageSource_h
#define __itkCOSMOSPointSpreadFunctionImageSource_h

#include "itkNumericTraits.h"
#include "itkParametricImageSource.h"

namespace itk
{

template< class TOutputImage >
class ITK_EXPORT COSMOSPointSpreadFunctionImageSource :
    public ParametricImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef COSMOSPointSpreadFunctionImageSource  Self;
  typedef ParametricImageSource< TOutputImage > Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /** Typedef for the output image PixelType. */
  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename OutputImageType::IndexType      OutputImageIndexType;
  typedef typename OutputImageType::RegionType     OutputImageRegionType;
  typedef typename OutputImageType::PointType      OutputImagePointType;
  typedef typename OutputImageType::PointValueType OutputImagePointValueType;
  typedef typename OutputImageType::SpacingType    OutputImageSpacingType;
  typedef typename OutputImageType::SizeType       OutputImageSizeType;
  typedef typename OutputImageType::SizeValueType  OutputImageSizeValueType;

  /** Parameter types. */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::ParametersValueType ParametersValueType;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(COSMOSPointSpreadFunctionImageSource, ParametricImageSource);

  /** Method for creation through the object factory. */
  //itkNewMacro(Self);

  /** Set/get  the magnification. */
  itkSetMacro(Magnification, double);
  itkGetConstMacro(Magnification, double);

  /** Set/get the numerical aperture (NA). */
  itkSetMacro(NumericalAperture, double);
  itkGetConstMacro(NumericalAperture, double);

  /** Set/get the emission wavelength. */
  itkSetMacro(EmissionWavelength, double);
  itkGetConstMacro(EmissionWavelength, double);

  /** Specify the design cover slip refractive index (unitless). */
  itkSetMacro(DesignCoverSlipRefractiveIndex, double);

  /** Get the design cover slip refractive index (unitless). */
  itkGetConstMacro(DesignCoverSlipRefractiveIndex, double);

  /** Specify the actual cover slip refractive index (unitless). */
  itkSetMacro(ActualCoverSlipRefractiveIndex, double);

  /** Get the actual cover slip refractive index (unitless). */
  itkGetConstMacro(ActualCoverSlipRefractiveIndex, double);

  /** Specify the design cover slip thickness (in micrometers). */
  itkSetMacro(DesignCoverSlipThickness, double);

  /** Get the design cover slip thickness (in micrometers). */
  itkGetConstMacro(DesignCoverSlipThickness, double);

  /** Specify the actual cover slip thickness (in micrometers). */
  itkSetMacro(ActualCoverSlipThickness, double);

  /** Get the actual cover slip thickness (in micrometers). */
  itkGetConstMacro(ActualCoverSlipThickness, double);

  /** Specify the design immersion oil refractive index (unitless). */
  itkSetMacro(DesignImmersionOilRefractiveIndex, double);

  /** Get the design immersion oil refractive index (unitless). */
  itkGetConstMacro(DesignImmersionOilRefractiveIndex, double);

  /** Specify the actual immersion oil refractive index (unitless). */
  itkSetMacro(ActualImmersionOilRefractiveIndex, double);

  /** Get the actual immersion oil refractive index (unitless). */
  itkGetConstMacro(ActualImmersionOilRefractiveIndex, double);

  /** Specify the design immersion oil thickness (in micrometers). */
  itkSetMacro(DesignImmersionOilThickness, double);

  /** Get the design immersion oil refractive index (in micrometers). */
  itkGetConstMacro(DesignImmersionOilThickness, double);

  /** Specify the design specimen layer refractive index (unitless). */
  itkSetMacro(DesignSpecimenLayerRefractiveIndex, double);

  /** Get the design specimen layer refractive index (unitless). */
  itkGetConstMacro(DesignSpecimenLayerRefractiveIndex, double);

  /** Specify the actual specimen layer refractive index (unitless). */
  itkSetMacro(ActualSpecimenLayerRefractiveIndex, double);

  /** Get the actual specimen layer refractive index (unitless). */
  itkGetConstMacro(ActualSpecimenLayerRefractiveIndex, double);

  /** Specify the actual point source depth in the specimen layer (in nanometers). */
  itkSetMacro(ActualPointSourceDepthInSpecimenLayer, double);

  /** Get the actual point source depth in the specimen layer (in nanometers). */
  itkGetConstMacro(ActualPointSourceDepthInSpecimenLayer, double);

  /** Get/set the shear in X. */
  itkSetMacro(ShearX, double);
  itkGetConstMacro(ShearX, double);

  /** Get/set the shear in Y. */
  itkSetMacro(ShearY, double);
  itkGetConstMacro(ShearY, double);

  /** Expects the parameters argument to contain values for ALL parameters. */
  virtual void SetParameters(const ParametersType& parameters);

  /** Gets the full parameters list. */
  virtual ParametersType GetParameters() const;

  /** Gets the total number of parameters. */
  virtual unsigned int GetNumberOfParameters() const;

protected:
  COSMOSPointSpreadFunctionImageSource();
  ~COSMOSPointSpreadFunctionImageSource();

  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  COSMOSPointSpreadFunctionImageSource(const COSMOSPointSpreadFunctionImageSource&); // purposely not implemented
  void operator=(const COSMOSPointSpreadFunctionImageSource&); // purposely not implemented

  // Unitless
  double m_Magnification;

  // Unitless
  double m_NumericalAperture;

  // Specified in nanometers
  double m_EmissionWavelength;

  // Unitless
  double m_DesignCoverSlipRefractiveIndex;
  double m_ActualCoverSlipRefractiveIndex;

  // Specified in micrometers
  double m_DesignCoverSlipThickness;
  double m_ActualCoverSlipThickness;

  // Unitless
  double m_DesignImmersionOilRefractiveIndex;
  double m_ActualImmersionOilRefractiveIndex;

  // Specified in micrometers
  double m_DesignImmersionOilThickness;

  // Unitless
  double m_DesignSpecimenLayerRefractiveIndex;
  double m_ActualSpecimenLayerRefractiveIndex;

  // Specified in nanometers
  double m_ActualPointSourceDepthInSpecimenLayer;

  // Specified in nanometers in X vs. nanometers in Z  
  double m_ShearX;

  // Specified in nanometers in Y vs. nanometers in Z
  double m_ShearY;

};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCOSMOSPointSpreadFunctionImageSource.hxx"
#endif

#endif

