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

#include "vtkWUWriter.h"
#include "wu/wuHeader.h"
#include <vtkErrorCode.h>
#include <vtkObjectFactory.h>
#include <vtkImageData.h>
#include <vtkPointData.h>

vtkCxxRevisionMacro(vtkWUWriter, "$Revision: 1.84 $")
vtkStandardNewMacro(vtkWUWriter)

//-----------------------------------------------------------------------------
// Constructor / Destructor
vtkWUWriter::vtkWUWriter()
{
}

vtkWUWriter::~vtkWUWriter()
{
}

//-----------------------------------------------------------------------------
// Print
void vtkWUWriter::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

void vtkWUWriter::WriteFileHeader(ofstream * file, vtkImageData * data) 
{
    int* extent = data->GetUpdateExtent();
    int dim[3] = {extent[1]-extent[0]+1,
                  extent[3]-extent[2]+1,
                  extent[5]-extent[4]+1};

    cosm::wuHeader header;
    header.sections( (unsigned int)dim[2] );
    header.rows( (unsigned int)dim[1] );
    header.columns( (unsigned int)dim[0] );
    
    switch (data->GetScalarType()) 
    {
        case VTK_DOUBLE: header.dataType( cosm::wuHeader::DOUBLE ); break;
        case VTK_FLOAT: header.dataType( cosm::wuHeader::FLOAT ); break;
        case VTK_INT: header.dataType( cosm::wuHeader::INT ); break;
        case VTK_SHORT: header.dataType( cosm::wuHeader::SHORT ); break;
        case VTK_UNSIGNED_SHORT: header.dataType( cosm::wuHeader::USHORT ); break;
        case VTK_UNSIGNED_CHAR: header.dataType( cosm::wuHeader::BYTE ); break;
        case VTK_CHAR: 
        case VTK_UNSIGNED_INT:
        case VTK_LONG:
        case VTK_UNSIGNED_LONG: break;
    }
    header.write(file);
    return;
}

void vtkWUWriter::WriteFile(ofstream *file, vtkImageData *data, int extent[6])
{
  int idx0, idx1, idx2;
  int rowLength; // in bytes
  void *ptr;
  unsigned long count = 0;
  unsigned long target;
  float progress = this->Progress;
  float area;
  int *wExtent;
  
  // Make sure we actually have data.
  if ( !data->GetPointData()->GetScalars())
    {
    vtkErrorMacro(<< "Could not get data from input.");
    return;
    }

  // take into consideration the scalar type
  switch (data->GetScalarType())
    {
    case VTK_UNSIGNED_CHAR:
      rowLength = sizeof(unsigned char); 
      break;
    case VTK_UNSIGNED_SHORT:
      rowLength = sizeof(unsigned short); 
      break;
    case VTK_SHORT:
      rowLength = sizeof(short); 
      break;
    case VTK_INT:
      rowLength = sizeof(int); 
      break;
    case VTK_FLOAT:
      rowLength = sizeof(float); 
      break;
    case VTK_DOUBLE:
      rowLength = sizeof(double); 
      break;
    default:
      vtkErrorMacro("WUWriter does not support type of scalars!");
      return; 
    }
  rowLength *= data->GetNumberOfScalarComponents();

  wExtent = this->GetInput()->GetWholeExtent();
  area = static_cast<float>(((extent[5] - extent[4] + 1)*(extent[3] - extent[2] + 1)*
                             (extent[1] - extent[0] + 1))) / 
         static_cast<float>(((wExtent[5] -wExtent[4] + 1)*(wExtent[3] -wExtent[2] + 1)*
                             (wExtent[1] -wExtent[0] + 1)));
    
  target = (unsigned long)((extent[5]-extent[4]+1)*
                           (extent[3]-extent[2]+1)/(50.0*area));
  target++;

  for (idx2 = extent[4]; idx2 <= extent[5]; ++idx2)
    {
    for (idx1 = extent[2]; idx1 <= extent[3]; idx1++)
      {
      if (!(count%target))
        {
        this->UpdateProgress(progress + count/(50.0*target));
        }
      count++;
      for (idx0 = extent[0]; idx0 <= extent[1]; idx0++)
        {
        ptr = data->GetScalarPointer(idx0, idx1, idx2);
        if ( ! file->write((char *)ptr, rowLength))
          {
          this->SetErrorCode(vtkErrorCode::OutOfDiskSpaceError);
          return;
          }
        }
      }
    }
}

