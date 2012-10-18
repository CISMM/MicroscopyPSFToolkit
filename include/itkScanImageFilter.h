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
#ifndef __itkScanImageFilter_h
#define __itkScanImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class ScanImageFilter
 * \brief Implements scan operations on an image along a selected direction.
 *
 * This class accumulates an image along a dimension, similar to
 * ProjectionImageFilter. Unlike ProjectionImageFilter, intermediate results of
 * the accumulation are stored in slices orthogonal to the scan direction. The
 * output image is the same size and dimensions.
 *
 * This class is parameterized over the type of the input and output images.
 *
 * Precision of the accumulation function is determined by the output type.
 *
 * \author Cory Quammen. Department of Computer Science, UNC Chapel Hill.
 */
template <class TInputImage, class TOutputImage, class TAccumulator>
class ITK_EXPORT ScanImageFilter :
  public ImageToImageFilter<TInputImage,TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ScanImageFilter                               Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScanImageFilter, ImageToImageFilter);

  /** Same convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::RegionType      InputImageRegionType;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef TOutputImage                             OutputImageType;
  typedef typename     OutputImageType::Pointer    OutputImagePointer;
  typedef typename     OutputImageType::RegionType OutputImageRegionType;
  typedef typename     OutputImageType::PixelType  OutputImagePixelType;

  typedef TAccumulator AccumulatorType;

  /** ImageDimension enumeration */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef enum {
    INCREASING_ORDER,
    DECREASING_ORDER
  } ScanOrder;

  /** Set/Get the direction in which to accumulate the data.  It must be set
   * before the update of the filter. Defaults to the last dimension. */
  itkSetMacro( ScanDimension, unsigned int );
  itkGetConstReferenceMacro( ScanDimension, unsigned int );

  /** Set/Get the order of traversal in the scan calculation. Increasing order
   * refers to increasing the index of the scan direction during the scan. */
  void SetScanOrderToIncreasing()
  {
    m_ScanOrder = INCREASING_ORDER;
    this->Modified();
  }

  void SetScanOrderToDecreasing()
  {
    m_ScanOrder = DECREASING_ORDER;
    this->Modified();
  }

protected:
  ScanImageFilter();
  virtual ~ScanImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Apply changes to the input image requested region. This method ensures
   * that the requested region extends to the beginning slice in the scan
   * direction, which depends on the scan order. */
  virtual void GenerateInputRequestedRegion();

  virtual void GenerateInputRequestedRegionForOutputRequestedRegion
    (const OutputImageRegionType& outputRegion,
     InputImageRegionType& inputRegion);

  /** Split the region into slices orthogonal to the scan direction. */
  virtual unsigned int SplitRequestedRegion(unsigned int i, unsigned int num, OutputImageRegionType& splitRegion);

  virtual void ThreadedGenerateData(
    const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId);

  virtual AccumulatorType NewAccumulator(unsigned long) const;

private:
  ScanImageFilter(const Self&); // purposely not implemented
  void operator=(const Self&); // purposely not implemented

  ScanOrder    m_ScanOrder;

  unsigned int m_ScanDimension;

}; // end class ScanImageFilter
} // end namespace itk

#include "itkScanImageFilter.hxx"

#endif // __itkScanImageFilter_h
