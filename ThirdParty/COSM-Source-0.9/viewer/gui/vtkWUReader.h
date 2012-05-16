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
                                                                                
#ifndef __vtkWUReader_h
#define __vtkWUReader_h

#include <vtkImageReader.h>
#include "wu/wuImage.h"
#include <string>

//-----------------------------------------------------------------------------
class VTK_EXPORT vtkWUReader : public vtkImageReader
{
public:
   static vtkWUReader *New();
   vtkTypeRevisionMacro(vtkWUReader, vtkImageReader);
   void PrintSelf(ostream& os, vtkIndent indent);
   virtual void SetFileName(const char *name);

  // Is the given file name a WU file?
  virtual int CanReadFile(const char* fname);

  // Get the file extensions for this format.
  // Returns a string with a space separated list of extensions in 
  // the format .extension
  virtual const char* GetFileExtensions()
  {
      return ".wu";
  }

  // Return a descriptive name for the file format that might be useful in a GUI.
  virtual const char* GetDescriptiveName()
  {
      return "Image with WU Header";
  }

protected:
   vtkWUReader();
   ~vtkWUReader();

   virtual void ExecuteInformation();
   virtual void ExecuteData(vtkDataObject *output);

   virtual void BuildData(vtkDataObject *output);
   virtual void LoadFileInformation();

private:

   void IncrementProgress(const unsigned long updateProgressTarget,
                          unsigned long &updateProgressCount);

   void LoadImageInMemory(unsigned char *dest,
                          const unsigned long updateProgressTarget,
                          unsigned long &updateProgressCount);

// Variables

   vtkTimeStamp fileTime;
   bool hasImage;

   //BTX
   // Number of columns of the image/volume to be loaded
   int NumColumns;
   // Number of lines of the image/volume to be loaded
   int NumLines;
   // Number of lines of the image/volume to be loaded
   int NumPlanes;
   // Total number of planes (or images) of the stack to be build.
   int TotalNumberOfPlanes;
   // Number of scalar components of the image to be loaded (1=monochrome 3=rgb)
   int NumComponents;
   // Type of the image[s]: 8/16/32 bits, signed/unsigned:
   std::string ImageType;
   // Pixel size (in number of bytes):
   size_t PixelSize;
   bool OwnFile;
   bool Execution;
   std::string ImageFileName;
   cosm::WUImage imageData;
  
   //ETX

};

//-----------------------------------------------------------------------------
#endif

