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
                                                                                
#ifndef VIEWER_BASE_H_
#define VIEWER_BASE_H_

#include "vtkImagePlaneWidget.h"
#include "vtkLineSource.h"
#include "vtkImageReader2.h"
#include <blitz/array.h>
#include <map>
#include <string>

class vtkWUReader;
class vtkTIFFReader;
class vtkPNGReader;
class vtkPNMReader;
class vtkJPEGReader;
class vtkRenderer;
class vtkImageMagnify;
class vtkImagePlaneWidget;
class vtkScalarBarActor;
class vtkActor;
class vtkRenderWindow;
class vtkWindowLevelLookupTable;

class ViewerBase
{

private:
  vtkImageReader2*	 Reader;
  vtkWUReader*	         WUReader;
  vtkTIFFReader*	 TIFFReader;
  vtkPNGReader*          PNGReader;
  vtkPNMReader*          PNMReader;
  vtkJPEGReader*         JPEGReader;
  vtkImagePlaneWidget**  Plane;
  vtkLineSource**        LineSource;
  vtkScalarBarActor*     ScalarBar;
  vtkActor**		 SliceActor;
  vtkActor**		 LineActor;
  vtkActor*		 Outline;
  vtkLookupTable* Table;
  vtkLookupTable* LogTable;
  vtkWindowLevelLookupTable* WindowLevelTable;
  std::string		 FileName;
  bool PipelineActive;
  bool RenderedPlane[3];
  bool RenderedSlice[3];
  bool RenderedVol;
  bool RenderedBar;

public:
  ViewerBase(void);

  ~ViewerBase()
    {
    }

  void AddActorsVol(
    vtkRenderWindow* aRenderWindow, 
    vtkRenderer* aRenderer
  );
  void AddActorsBar(
    vtkRenderer* aRenderer
  );
  void AddActorsPlane(
    int index, 
    vtkRenderWindow* aRenderWindow, 
    vtkRenderer* aRenderer
  );
  void AddActorsSlice(
    int index, 
    vtkRenderWindow* aRenderWindow, 
    vtkRenderer* aRenderer
  );

  // Set/Get the current slice to display

  int GetZSlice() {return GetSlice(2); };
  void SetZSlice(int s) { SetSlice(2,s); };
  int GetYSlice() {return GetSlice(1); };
  void SetYSlice(int s) { SetSlice(1,s); };
  int GetXSlice() {return GetSlice(0); };
  void SetXSlice(int s) { SetSlice(0,s); };

  int GetZSize( void ) 
  { 
    if ( Reader == NULL ) return 0;
    int* extent = Reader->GetDataExtent(); 
    return extent[5];
  }
  int GetYSize( void )
  { 
    if ( Reader == NULL ) return 0;
    int* extent = Reader->GetDataExtent(); 
    return extent[3];
  }
  int GetXSize( void )
  { 
    if ( Reader == NULL ) return 0;
    int* extent = Reader->GetDataExtent(); 
    return extent[1];
  }

  void GetLine( double* line, int z, int y, int x, int sel );

  void SetLevel(
      double val
  );

  void SetWindow(
      double val
  );
  void SetScale(
      int val
  );

  void SetFileName( const std::string& filename ) { this->FileName = filename; };
  void WriteFile( const std::string& filename );
  void Pipeline(void);

protected:
  int GetSlice(int i);
  void SetSlice(int i, int s);
  bool CheckReader( void );

private:
  ViewerBase (const ViewerBase&); // Not Implemented.
  ViewerBase& operator= (const ViewerBase&); // Not Implemented.

};

#endif /* VIEWER_BASE_H_ */
