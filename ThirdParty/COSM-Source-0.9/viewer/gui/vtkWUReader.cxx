/****************************************************************************
 * Copyright (c) 2007 Einir Valdimarsson and Chrysanthe Preza
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 ****************************************************************************/
                                                                                
#include "vtkWUReader.h"
#include "wu/wuHeader.h"
#include <vtkObjectFactory.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <blitz/array.h>

vtkCxxRevisionMacro(vtkWUReader, "$Revision: 1.84 $")
vtkStandardNewMacro(vtkWUReader)

//-----------------------------------------------------------------------------
// Constructor / Destructor
vtkWUReader::vtkWUReader()
{
   this->OwnFile = true;
   this->Execution = false;
}

vtkWUReader::~vtkWUReader()
{
}

//-----------------------------------------------------------------------------
// Print
void vtkWUReader::PrintSelf(ostream &os, vtkIndent indent)
{
   this->Superclass::PrintSelf(os,indent);
}

//-----------------------------------------------------------------------------
// Public

/*
 * Sets up a filename to be read.
 */
void vtkWUReader::SetFileName(const char *name) 
{
   vtkImageReader2::SetFileName(name);
   ImageFileName = name;
   this->Modified();
}

//-----------------------------------------------------------------------------
// Protected
/*
 * Configure the output e.g. WholeExtent, spacing, origin, scalar type...
 */
void vtkWUReader::ExecuteInformation()
{
   if (this->Execution)
      return;

   this->Execution=true;

   if (this->MTime > this->fileTime)
   {
      this->TotalNumberOfPlanes = 0;
      this->LoadFileInformation();

      if ( this->TotalNumberOfPlanes == 0)
      {
         //vtkErrorMacro(<< "File set is not valid. Exiting...");
         return;
      }

      // if the user has not set the extent, but has set the VOI
      // set the z axis extent to the VOI z axis
      if (this->DataExtent[4]==0 && this->DataExtent[5] == 0 &&
         (this->DataVOI[4] || this->DataVOI[5]))
      {
         this->DataExtent[4] = this->DataVOI[4];
         this->DataExtent[5] = this->DataVOI[5];
      }

      // When the user has set the VOI, check it's coherence with the file content.
      if (this->DataVOI[0] || this->DataVOI[1] || 
      this->DataVOI[2] || this->DataVOI[3] ||
      this->DataVOI[4] || this->DataVOI[5])
      { 
         if ((this->DataVOI[0] < 0) ||
             (this->DataVOI[1] >= this->NumColumns) ||
             (this->DataVOI[2] < 0) ||
             (this->DataVOI[3] >= this->NumLines) ||
             (this->DataVOI[4] < 0) ||
             (this->DataVOI[5] >= this->TotalNumberOfPlanes ))
         {
            vtkWarningMacro(<< "The requested VOI is larger than expected extent.");
            this->DataVOI[0] = 0;
            this->DataVOI[1] = this->NumColumns > 0 ? this->NumColumns - 1 : 0;
            this->DataVOI[2] = 0;
            this->DataVOI[3] = this->NumLines > 0 ? this->NumLines - 1 : 0;
            this->DataVOI[4] = 0;
            this->DataVOI[5] = this->TotalNumberOfPlanes > 0 ? this->TotalNumberOfPlanes - 1 : 0;
         }
      }

      // Set the Extents.
      this->DataExtent[0] = 0;
      this->DataExtent[1] = this->NumColumns > 0 ? this->NumColumns - 1 : 0;
      this->DataExtent[2] = 0;
      this->DataExtent[3] = this->NumLines > 0 ? this->NumLines - 1 : 0;
      this->DataExtent[4] = 0;
      this->DataExtent[5] = this->TotalNumberOfPlanes > 0 ? this->TotalNumberOfPlanes - 1 : 0;
  
      // We do need to set the Endian related stuff (by using
      // this->SetDataByteOrderToBigEndian() or SetDataByteOrderToLittleEndian()
      // TODO - set Endian 

      // We do need to set up the data type for downstream filters:
      if      ( ImageType == "8U" )
      {
         vtkDebugMacro(<< "8 bits unsigned image");
         this->SetDataScalarTypeToUnsignedChar(); 
      }
      else if ( ImageType == "8S" )
      {
         vtkErrorMacro(<< "Cannot handle 8 bit signed files");
         return;
      }
      else if ( ImageType == "16U" )
      {
         vtkDebugMacro(<< "16 bits unsigned image");
         this->SetDataScalarTypeToUnsignedShort();
      }
      else if ( ImageType == "16S" )
      {
         vtkDebugMacro(<< "16 bits signed image");
         this->SetDataScalarTypeToShort();
      }
      else if ( ImageType == "32U" )
      {
         vtkDebugMacro(<< "32 bits unsigned image");
         vtkDebugMacro(<< "WARNING: forced to signed int !");
         this->SetDataScalarTypeToInt();
      }
      else if ( ImageType == "32S" )
      {
         vtkDebugMacro(<< "32 bits signed image");
         this->SetDataScalarTypeToInt();
      }
      else if ( ImageType == "32F" )
      {
         vtkDebugMacro(<< "32 bits float image");
         this->SetDataScalarTypeToFloat();
      }
      else if ( ImageType == "FD" )
      {
         vtkDebugMacro(<< "64 bits Double image");
         this->SetDataScalarTypeToDouble();
      }
      //Set number of scalar components:
      this->SetNumberOfScalarComponents(this->NumComponents);

      this->fileTime = this->MTime;
   }

   this->Superclass::ExecuteInformation();

   this->GetOutput()->SetUpdateExtentToWholeExtent();

   this->Execution=false;
}
 
/*
 * Update => ouput->Update => UpdateData => Execute => ExecuteData 
 * (see vtkSource.cxx for last step).
 * This function (redefinition of vtkImageReader::ExecuteData, see 
 * VTK/IO/vtkImageReader.cxx) reads a data from a file. The data
 * extent/axes are assumed to be the same as the file extent/order.
 */
void vtkWUReader::ExecuteData(vtkDataObject *output)
{
   vtkImageReader::ExecuteData(output);
   vtkImageData *data=vtkImageData::SafeDownCast(output);
   data->SetExtent(this->DataExtent);
   this->BuildData(this->GetOutput());
}

void vtkWUReader::BuildData(vtkDataObject *output)
{
   vtkImageData *data = this->AllocateOutputData(output);
   data->GetPointData()->GetScalars()->SetName("wuImage-Volume");

   // Test if output has valid extent
   // Prevent memory errors
   if((this->DataExtent[1]-this->DataExtent[0]>=0) &&
      (this->DataExtent[3]-this->DataExtent[2]>=0) &&
      (this->DataExtent[5]-this->DataExtent[4]>=0))
   {
      // Variables for the UpdateProgress. We shall use 50 steps to signify
      // the advance of the process:
      unsigned long UpdateProgressTarget = (unsigned long) ceil (this->NumLines
                                         * this->TotalNumberOfPlanes
                                         / 50.0);
      // The actual advance measure:
      unsigned long UpdateProgressCount = 0;

      // Filling the allocated memory space with each image/volume:
      unsigned char *Dest = (unsigned char *)data->GetScalarPointer();
      this->LoadImageInMemory(Dest, UpdateProgressTarget, UpdateProgressCount); 
   }
}

/*
 * Read the image 
 */
void vtkWUReader::LoadFileInformation()
{
   this->imageData.ReadData( this->ImageFileName );

   this->OwnFile=true;

   // Get the image caracteristics
   if ( this->imageData.IsFloat() )
   {
       vtkDebugMacro(<< "Data type: float");
       blitz::Array<float,3> image = this->imageData.GetFloatArray();
       this->NumColumns = image.extent(2);
       this->NumLines   = image.extent(1);
       this->NumPlanes  = image.extent(0);
       this->ImageType = "32F";
       this->PixelSize = 4;
   }
   else if ( this->imageData.IsDouble() )
   {
       vtkDebugMacro(<< "Data type: double");
       blitz::Array<double,3> image = this->imageData.GetDoubleArray();
       this->NumColumns = image.extent(2);
       this->NumLines   = image.extent(1);
       this->NumPlanes  = image.extent(0);
       this->ImageType = "FD";
       this->PixelSize = 8;
   }
   else if ( this->imageData.IsLongDouble() )
   {
       vtkDebugMacro(<< "Data type: long double");
       this->imageData.ConvertToDouble();
       // no long double image format - convert to double for display
       blitz::Array<double,3> image = this->imageData.GetDoubleArray();
       this->NumColumns = image.extent(2);
       this->NumLines   = image.extent(1);
       this->NumPlanes  = image.extent(0);
       this->ImageType = "FD";
       this->PixelSize = 8;
   }
   else if ( this->imageData.IsUshort() )
   {
       vtkDebugMacro(<< "Data type: unsigned short");
       blitz::Array<unsigned short,3> image = this->imageData.GetUshortArray();
       this->NumColumns = image.extent(2);
       this->NumLines   = image.extent(1);
       this->NumPlanes  = image.extent(0);
       this->ImageType = "16U";
       this->PixelSize = 2;
   }

   this->TotalNumberOfPlanes = this->NumPlanes;


   this->DataOrigin[0] = 0;
   this->DataOrigin[1] = 0;
   this->DataOrigin[2] = 0;

   this->DataSpacing[0] = 1;
   this->DataSpacing[1] = 1;
   this->DataSpacing[2] = 1;

   // Set the image data caracteristics
   this->NumComponents = 1;
}

//-----------------------------------------------------------------------------
// Private

void vtkWUReader::IncrementProgress(const unsigned long updateProgressTarget,
                                      unsigned long &updateProgressCount)
{
   // Update progress related for bad files:
   updateProgressCount += this->NumLines;
   if (updateProgressTarget > 0)
   {
      if (!(updateProgressCount%updateProgressTarget))
      {
         this->UpdateProgress(updateProgressCount/(50.0*updateProgressTarget));
      }
   }
}

/*
 * Loads the contents of the image/volume atthe Dest memory address. 
 * Returns the size of the data loaded.
 * \remarks Assume that if (f != NULL) then its caracteristics match
 * with the previous ones
 */
void vtkWUReader::LoadImageInMemory(
             unsigned char *dest,
             const unsigned long updateProgressTarget,
             unsigned long &updateProgressCount)
{
   vtkDebugMacro(<< "NumComponents:" << NumComponents);
   vtkDebugMacro(<< "Copying to memory image [" << this->ImageFileName.c_str() << "]");

   // If the data structure of vtk for image/volume representation
   // were straigthforwards the following would be enough.
   // But vtk chooses to invert the lines of an image, that is the last
   // line comes first (for some axis related reasons?). Hence we need
   // to load the image line by line, starting from the end.

   unsigned char *src = 0;
   int lineSize   = this->NumComponents * this->NumColumns * this->PixelSize;
   int planeSize  = lineSize * this->NumLines;
   if ( this->imageData.IsFloat() )
   {
       src = (unsigned char*)(this->imageData.GetFloatData());  
   } 
   else if ( this->imageData.IsDouble() )
   {
       src = (unsigned char*)(this->imageData.GetDoubleData());  
   }
   else if ( this->imageData.IsUshort() )
   {
       src = (unsigned char*)(this->imageData.GetUshortData());  
   }
   else
   {
       vtkDebugMacro(<< "Image data type unknown");
       return;
   }

//   unsigned char *dst = dest + planeSize - lineSize;
   unsigned char *dst = dest;
   for (int plane = 0; plane < this->NumPlanes; plane++)
   {
      for (int line = 0; line < this->NumLines; line++)
      {
         // Copy one line at proper destination:
         memcpy((void*)dst, (void*)src, lineSize);
         src += lineSize;
         dst += lineSize;
         //dst -= lineSize;
         // Update progress related:
         if (!(updateProgressCount%updateProgressTarget))
         {
            this->UpdateProgress(updateProgressCount/(50.0*updateProgressTarget));
         }
         updateProgressCount++;
      }
      //dst += 2 * planeSize;
   }
}

int vtkWUReader::CanReadFile(const char* fname) 
{
    cosm::wuHeader header;
    int res = header.read(fname);
    return  res == -1 ? 0 : 3;
}

//-----------------------------------------------------------------------------
