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

#ifndef __itkHaeberleCOSMOSPointSpreadFunctionImageSource_h
#define __itkHaeberleCOSMOSPointSpreadFunctionImageSource_h

#include "itkCOSMOSPointSpreadFunctionImageSource.h"
#include "itkNumericTraits.h"

#include <complex>

/* Unfortunately, this line is required to use the COSMOS PSF source. */
using std::complex;

#include "psf/functor.h"
#include "psf/haeberlePsfFunctor.h"

namespace itk
{

/** \class HaeberleCOSMOSPointSpreadFunctionImageSource
 * \brief Generate a synthetic point-spread function according to the
 * Haeberle model.
 *
 * The Haeberle point-spread function model is based on the vectorial
 * model of light propagation in widefield fluorescence
 * microscopes. This image source generates images according to this
 * IMPORTANT: Please pay attention to the units each method
 * expects. Some take nanometers, some take micrometers, and some
 * take millimeters.
 *
 * \ingroup DataSources Multithreaded
 */
template <class TOutputImage>
class ITK_EXPORT HaeberleCOSMOSPointSpreadFunctionImageSource :
    public COSMOSPointSpreadFunctionImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef HaeberleCOSMOSPointSpreadFunctionImageSource         Self;
  typedef COSMOSPointSpreadFunctionImageSource< TOutputImage > Superclass;
  typedef SmartPointer< Self >                                 Pointer;
  typedef SmartPointer< const Self >                           ConstPointer;

  /** Typedef for the output image PixelType. */
  typedef TOutputImage                                   OutputImageType;
  typedef typename Superclass::OutputImagePixelType      OutputImagePixelType;
  typedef typename Superclass::OutputImageIndexType      OutputImageIndexType;
  typedef typename Superclass::OutputImageRegionType     OutputImageRegionType;
  typedef typename Superclass::OutputImagePointType      OutputImagePointType;
  typedef typename Superclass::OutputImagePointValueType OutputImagePointValueType;
  typedef typename Superclass::OutputImageSpacingType    OutputImageSpacingType;
  typedef typename Superclass::OutputImageSizeType       OutputImageSizeType;
  typedef typename Superclass::OutputImageSizeValueType  OutputImageSizeValueType;

  /** Parameter types. */
  typedef typename Superclass::ParametersType      ParametersType;
  typedef typename Superclass::ParametersValueType ParametersValueType;

  itkStaticConstMacro(ImageDimension, unsigned int,
		      TOutputImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(HaeberleCOSMOSPointSpreadFunctionImageSource, ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

protected:
  HaeberleCOSMOSPointSpreadFunctionImageSource();
  ~HaeberleCOSMOSPointSpreadFunctionImageSource();

  /** The COSM Haeberle computations are not thread-safe, so we
   * override the single-threaded GenerateData() here. */
  void GenerateData();

  /** Computes the light intensity at a specified point. */
  double ComputeSampleValue(const OutputImagePointType & point);

private:
  HaeberleCOSMOSPointSpreadFunctionImageSource(const HaeberleCOSMOSPointSpreadFunctionImageSource&); //purposely not implemented
  void operator=(const HaeberleCOSMOSPointSpreadFunctionImageSource&); //purposely not implemented

  /** The Haeberle functor that does all the work. */
  cosm::HaeberlePsfFunctor< double > * m_HaeberleFunctor;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHaeberleCOSMOSPointSpreadFunctionImageSource.hxx"
#endif

#endif
