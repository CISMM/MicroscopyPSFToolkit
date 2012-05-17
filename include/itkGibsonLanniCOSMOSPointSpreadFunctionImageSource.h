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

#ifndef __itkGibsonLanniCOSMOSPointSpreadFunctionImageSource_h
#define __itkGibsonLanniCOSMOSPointSpreadFunctionImageSource_h

#include "itkCOSMOSPointSpreadFunctionImageSource.h"
#include "itkNumericTraits.h"

#include <complex>

/* Unfortunately, this line is required to use the COSMOS PSF source. */
using std::complex;

#include "psf/functor.h"
#include "psf/gibsonLaniPsfFunctor.h"

namespace itk
{

/** \class GibsonLanniCOSMOSPointSpreadFunctionImageSource
 * \brief Generate a synthetic point-spread function according to the
 * Gibson-Lanni model.
 *
 * The Gibson-Lanni point-spread function model takes into account optical
 * path differences from the design conditions of an objective in a
 * widefield fluorescence microscope. This image source generates images
 * according to this model. IMPORTANT: Please pay attention to the units
 * each method expects. Some take nanometers, some take micrometers, and some
 * take millimeters.
 *
 * \ingroup DataSources Multithreaded
 */
template< class TOutputImage >
class ITK_EXPORT GibsonLanniCOSMOSPointSpreadFunctionImageSource :
    public COSMOSPointSpreadFunctionImageSource< TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef GibsonLanniCOSMOSPointSpreadFunctionImageSource      Self;
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
  itkTypeMacro(GibsonLanniCOSMOSPointSpreadFunctionImageSource, ParametricImageSource);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

protected:
  GibsonLanniCOSMOSPointSpreadFunctionImageSource();
  ~GibsonLanniCOSMOSPointSpreadFunctionImageSource();

  /** I made changes to the integrators in cquadpack so that they
   *  could be safely used by multiple threads. */
  void ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread,
                            ThreadIdType threadId);

private:
  GibsonLanniCOSMOSPointSpreadFunctionImageSource(const GibsonLanniCOSMOSPointSpreadFunctionImageSource&); //purposely not implemented
  void operator=(const GibsonLanniCOSMOSPointSpreadFunctionImageSource&); //purposely not implemented
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGibsonLanniCOSMOSPointSpreadFunctionImageSource.hxx"
#endif

#endif
