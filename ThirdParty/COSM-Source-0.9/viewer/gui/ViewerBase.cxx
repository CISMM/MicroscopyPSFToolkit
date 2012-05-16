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
                                                                                
#include "ViewerBase.h"
// VTK Graphics
#include "vtkWUReader.h"
#include "vtkTIFFReader.h"
#include "vtkPNGReader.h"
#include "vtkPNMReader.h"
#include "vtkJPEGReader.h"
#include "vtkWUWriter.h"
#include <vtkPlaneSource.h>
#include <vtkCamera.h>
#include <vtkLightKit.h>
#include <vtkImageData.h>

// VTK Rendering
#include "vtkWindowLevelLookupTable.h"
#include "vtkImageMapper.h"
#include "vtkActor2D.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkImageReslice.h"
#include "vtkImagePlaneWidget.h"
#include "vtkScalarBarActor.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkCellPicker.h"
#include "vtkProperty.h"

// ----------------------------------------------------------------------------
//      C o s m    V i e w e r
// ----------------------------------------------------------------------------

ViewerBase::ViewerBase(void) 
: Reader(0), WUReader(0), TIFFReader(0), PNGReader(0), PNMReader(0), JPEGReader(0), Plane(0), LineSource(0), ScalarBar(0), SliceActor(0), LineActor(0), Outline(0), Table(0), LogTable(0), WindowLevelTable(0), PipelineActive(false), RenderedVol(false), RenderedBar(false)
{
    for ( int i = 0; i < 3; i++ )
    {
        RenderedPlane[i] = false;
        RenderedSlice[i] = false;
    }
    WUReader = vtkWUReader::New();
    TIFFReader = vtkTIFFReader::New();
    PNGReader = vtkPNGReader::New();
    PNMReader = vtkPNMReader::New();
    JPEGReader = vtkJPEGReader::New();
}

void
ViewerBase::Pipeline(void)
{
    if ( !CheckReader() ) return;
    Reader->Update();
    if ( PipelineActive )
    {
        vtkFloatingPointType* range = Reader->GetOutput()->GetScalarRange();
        WindowLevelTable->SetLevel(0.5*(range[1] + range[0]));
        WindowLevelTable->SetWindow(1.0*(range[1] - range[0]));
        for ( int i = 0; i < 3; i++ )
        {
            if ( Plane[i] != NULL )
            {
	        Plane[i]->SetInput((vtkDataSet*)Reader->GetOutput());
            }
            vtkOutlineFilter* outlineFilter = vtkOutlineFilter::New();
            outlineFilter->SetInput((vtkDataSet*)Reader->GetOutput());
            vtkPolyDataMapper* outlineMapper = vtkPolyDataMapper::New();
            outlineMapper->SetInput(outlineFilter->GetOutput());
            if ( Outline != NULL )
            {
                Outline->SetMapper(outlineMapper);
            }
            outlineFilter->Delete();
            outlineMapper->Delete();
        }
	return;
    }
    else if ( Reader == NULL )
    {
        return;
    }

    {
        vtkFloatingPointType* range = Reader->GetOutput()->GetScalarRange();

        WindowLevelTable = vtkWindowLevelLookupTable::New();
        WindowLevelTable->SetLevel(0.5*(range[1] + range[0]));
        WindowLevelTable->SetWindow(1.0*(range[1] - range[0]));
        WindowLevelTable->Build();

        LogTable = vtkLookupTable::New();
        LogTable->SetScaleToLog10();

        Table = WindowLevelTable;

        ScalarBar = vtkScalarBarActor::New();
        ScalarBar->SetLookupTable(Table);
        ScalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
        ScalarBar->GetPositionCoordinate()->SetValue(0.1,0.02);
        ScalarBar->SetWidth(0.9);
        ScalarBar->SetHeight(1.12);

        vtkCellPicker* picker = vtkCellPicker::New();
	picker->SetTolerance(0.005);

        Plane = new vtkImagePlaneWidget*[3];
        for ( int i = 0; i < 3; i++ )
        {
            Plane[i] = vtkImagePlaneWidget::New();
	    Plane[i]->DisplayTextOn();
	    Plane[i]->SetInput((vtkDataSet*)Reader->GetOutput());
	    Plane[i]->SetSliceIndex(0);
	    Plane[i]->SetPicker(picker);
	    Plane[i]->GetPlaneProperty()->SetColor(1,0,0);
            Plane[i]->SetLookupTable(Table);
        }
        Plane[2]->SetPlaneOrientationToZAxes();
        Plane[1]->SetPlaneOrientationToYAxes();
        Plane[0]->SetPlaneOrientationToXAxes();

        vtkOutlineFilter* outlineFilter = vtkOutlineFilter::New();
        outlineFilter->SetInput((vtkDataSet*)Reader->GetOutput());
        vtkPolyDataMapper* outlineMapper = vtkPolyDataMapper::New();
        outlineMapper->SetInput(outlineFilter->GetOutput());
        Outline = vtkActor::New();
        Outline->SetMapper(outlineMapper);

        SliceActor = new vtkActor*[3];

        for ( int i = 0; i < 3; i++ )
        {
            vtkPlaneSource* slicePlane = vtkPlaneSource::New();
            if ( i == 1 )
            {
                slicePlane->SetOrigin(0.0,0.0,0.0);
                slicePlane->SetPoint1(0.0,1.0,0.0);
                slicePlane->SetPoint2(1.0,0.0,0.0);
            }
            else
            {
                slicePlane->SetOrigin(0.0,0.0,0.0);
                slicePlane->SetPoint1(1.0,0.0,0.0);
                slicePlane->SetPoint2(0.0,1.0,0.0);
            }
            vtkPolyDataMapper* sliceMapper = vtkPolyDataMapper::New();
            sliceMapper->SetInput(slicePlane->GetOutput());
            sliceMapper->ImmediateModeRenderingOn();

            SliceActor[i] = vtkActor::New();
            SliceActor[i]->SetTexture(Plane[i]->GetTexture());
            SliceActor[i]->SetMapper(sliceMapper);

            slicePlane->Delete();
            sliceMapper->Delete();

        }

        LineActor = new vtkActor*[6];
        LineSource = new vtkLineSource*[6];
        for ( int i = 0; i < 6; i++ )
        {
            LineSource[i] = vtkLineSource::New();
            if ( i < 3 && i == 1 )
            {
                LineSource[i]->SetPoint1(0.5, 0.0, 0.0);
                LineSource[i]->SetPoint2(0.5, 1.0, 0.0);
            }
            else if ( i < 3 && i != 1 )
            {
                LineSource[i]->SetPoint1(0.0, 0.5, 0.0);
                LineSource[i]->SetPoint2(1.0, 0.5, 0.0);
            }
            else if ( i >= 3 && i == 1 )
            {
                LineSource[i]->SetPoint1(0.0, 0.5, 0.0);
                LineSource[i]->SetPoint2(1.0, 0.5, 0.0);
            }
            else if ( i >= 3 && i != 1 )
            {
                LineSource[i]->SetPoint1(0.5, 0.0, 0.0);
                LineSource[i]->SetPoint2(0.5, 1.0, 0.0);
            }
            vtkPolyDataMapper* lineMapper = vtkPolyDataMapper::New();
            lineMapper->SetInput(LineSource[i]->GetOutput() );
            LineActor[i] = vtkActor::New();
            LineActor[i]->SetMapper(lineMapper); 
            LineActor[i]->GetProperty()->SetOpacity(0.5);
           
            lineMapper->Delete();
        }

	PipelineActive = true;
        picker->Delete();
        outlineFilter->Delete();
        outlineMapper->Delete();
    }
}

void  
ViewerBase::AddActorsVol (
    vtkRenderWindow* aRenderWindow, 
    vtkRenderer* aRenderer
) {
    if ( !PipelineActive )
    {
	return;
    }
    if ( aRenderer == NULL || Plane == NULL )
    {
	return;
    }
    if ( RenderedVol )
    {
        for ( int i = 0; i < 3; i++ )
        {
            if ( Plane[i] != NULL )
            {
               vtkRenderWindowInteractor* vtkInteractor = aRenderWindow->GetInteractor();
                vtkInteractorStyleTrackballCamera* vtkTrackball = vtkInteractorStyleTrackballCamera::New();
                vtkInteractor->SetInteractorStyle(vtkTrackball);
                vtkTrackball->Delete();
               Plane[i]->SetInteractor(vtkInteractor);
               Plane[i]->On();
               Plane[i]->UpdatePlacement();
               aRenderer->ResetCamera();
            }
        }
        return;
    }
    if ( Outline != NULL )
    {
        aRenderer->AddActor(Outline);
        aRenderer->ResetCamera();
        aRenderer->SetBackground(0,0,0.5);
        RenderedVol = true;
    }
    for ( int i = 0; i < 3; i++ )
    {
        if ( Plane[i] != NULL )
        {
            Plane[i]->SetInteractor(aRenderWindow->GetInteractor());
            Plane[i]->On();
        }
    }
}

void  
ViewerBase::AddActorsBar(
    vtkRenderer* aRenderer
) {
    if ( !PipelineActive || RenderedBar )
    {
	return;
    }
    if ( aRenderer == NULL )
    {
	return;
    }
    if ( ScalarBar != NULL )
    {
	aRenderer->AddActor(ScalarBar); 
        aRenderer->SetBackground(0,0,0.5);
        aRenderer->ResetCamera();
        RenderedBar = true;
    }
}

void  
ViewerBase::AddActorsPlane(
    int index, 
    vtkRenderWindow* aRenderWindow, 
    vtkRenderer* aRenderer
) {
    if ( !PipelineActive || RenderedPlane[index] )
    {
	return;
    }
    if ( aRenderer == NULL )
    {
	return;
    }
    if ( SliceActor != NULL && SliceActor[index] != NULL )
    {
        aRenderer->AddActor(SliceActor[index]);
        aRenderer->AddActor(LineActor[index]);
        aRenderer->AddActor(LineActor[index+3]);
        aRenderer->SetBackground(0,0,0.5);
        aRenderer->ResetCamera();
	RenderedPlane[index] = true;
    }
}

void  
ViewerBase::AddActorsSlice(
    int index, 
    vtkRenderWindow* aRenderWindow, 
    vtkRenderer* aRenderer
) {
    if ( !PipelineActive || RenderedSlice[index] )
    {
	return;
    }
    if ( aRenderer == NULL )
    {
	return;
    }
    if ( SliceActor != NULL && SliceActor[index] != NULL )
    {
        aRenderer->AddActor(SliceActor[index]);
        aRenderer->SetBackground(0,0,0.5);
        aRenderer->ResetCamera();
	RenderedSlice[index] = true;
    }
}

void ViewerBase::SetLevel(
    double val
) {
    if ( WindowLevelTable )
    {
        vtkFloatingPointType* range = Reader->GetOutput()->GetScalarRange();
        this->WindowLevelTable->SetLevel (val * (range[1] + range[0]));
    }
}
                                                                                
void ViewerBase::SetWindow(
    double val
) {
    if ( WindowLevelTable )
    { 
        vtkFloatingPointType* range = Reader->GetOutput()->GetScalarRange();
        this->WindowLevelTable->SetWindow (val * (range[1] - range[0]));
    }
}

void ViewerBase::SetScale(
    int val
) {
    if ( Table )
    {
        if ( val == 0 ) 
        { 
            Table = WindowLevelTable;
        }
        else
        {
            Table = LogTable;
        }
        ScalarBar->SetLookupTable(Table);
        for ( int i = 0; i < 3; i++ )
        {
            Plane[i]->SetLookupTable(Table);
        }
    }
}

void ViewerBase::SetSlice(
    int i, 
    int s
) {
    if (Plane != 0 &&  Plane[i] != 0 )
    {
        this->Plane[i]->SetSliceIndex(s);
    }
    if ( LineSource != 0 )
    {
        if ( i == 2 && LineSource[1] != 0 && LineSource[0] != 0 )
        {
            double val = double(s)/GetZSize();
            LineSource[0]->SetPoint1(0.0, val, 0.0); 
            LineSource[0]->SetPoint2(1.0, val, 0.0); 
            LineSource[4]->SetPoint1(0.0, val, 0.0); 
            LineSource[4]->SetPoint2(1.0, val, 0.0); 
            LineActor[0]->GetProperty()->SetColor(0.0, 1.0, 1.0);
            LineActor[4]->GetProperty()->SetColor(1.0, 0.0, 0.0);
        }
        if ( i == 1 && LineSource[2] != 0 && LineSource[3] != 0 )
        {
            double val = double(s)/GetYSize();
            LineSource[2]->SetPoint1(0.0, val, 0.0); 
            LineSource[2]->SetPoint2(1.0, val, 0.0); 
            LineSource[3]->SetPoint1(val, 0.0, 0.0); 
            LineSource[3]->SetPoint2(val, 1.0, 0.0); 
            LineActor[2]->GetProperty()->SetColor(1.0, 0.0, 0.0);
            LineActor[3]->GetProperty()->SetColor(0.0, 1.0, 0.0);
        }
        if ( i == 0 && LineSource[4] != 0 && LineSource[5] != 0 )
        {
            double val = double(s)/GetXSize();
            LineSource[1]->SetPoint1(val, 0.0, 0.0); 
            LineSource[1]->SetPoint2(val, 1.0, 0.0); 
            LineSource[5]->SetPoint1(val, 0.0, 0.0); 
            LineSource[5]->SetPoint2(val, 1.0, 0.0); 
            LineActor[1]->GetProperty()->SetColor(0.0, 1.0, 0.0);
            LineActor[5]->GetProperty()->SetColor(0.0, 1.0, 1.0);
        }
    }
}

int ViewerBase::GetSlice(
    int i
) {
    return (Plane != 0 && Plane[i] != 0) ? this->Plane[i]->GetSliceIndex() : 0;
}


void ViewerBase::GetLine( double* profile, int z, int y, int x, int val )
{
    if ( Reader != NULL )
    {
        int* extent = Reader->GetDataExtent();

        int zMin = (val != 2) ? z : extent[4];
        int zMax = (val != 2) ? z : extent[5];
        int yMin = (val != 1) ? y : extent[2];
        int yMax = (val != 1) ? y : extent[3];
        int xMin = (val != 0) ? x : extent[0];
        int xMax = (val != 0) ? x : extent[1];

        int dim[3] = {extent[1]-extent[0]+1,
                      extent[3]-extent[2]+1,
                      extent[5]-extent[4]+1};

        vtkImageData *data=vtkImageData::SafeDownCast(Reader->GetOutput());

        int columnSize = data->GetNumberOfScalarComponents() * data->GetScalarSize();
        int lineSize = columnSize * dim[0];
        int planeSize =  lineSize * dim[1];

        unsigned char *src = (unsigned char *)data->GetScalarPointerForExtent(extent);
        int i = 0;
        for (int plane = zMin; plane <= zMax; plane++)
        {
            for (int line = yMin; line <= yMax; line++)
            {
                for (int column = xMin; column <= xMax; column++)
                {
                    int offset = plane*planeSize + line*lineSize + column*columnSize;
                    unsigned char* ptr = src + offset;
                    switch ( data->GetScalarSize() ) 
                    {
                        case  2: profile[i++] = (double)(*((unsigned short*)(ptr))); break;
                        case  4: profile[i++] = (double)(*((float*)(ptr))); break;
                        case  8: profile[i++] = (double)(*((double*)(ptr))); break;
                        case 16: profile[i++] = (double)(*((long double*)(ptr))); break;
                    }
                }
            }
        }
    }
}

bool ViewerBase::CheckReader( void ) {
   if ( WUReader->CanReadFile(this->FileName.c_str()) ) 
   {
       Reader = WUReader;
   } 
   else if ( TIFFReader->CanReadFile(this->FileName.c_str()) ) 
   {
       Reader = TIFFReader; 
   }
   else if ( PNGReader->CanReadFile(this->FileName.c_str()) ) 
   {
       Reader = PNGReader; 
   }
   else if ( PNMReader->CanReadFile(this->FileName.c_str()) ) 
   {
       Reader = PNMReader; 
   }
   else if ( JPEGReader->CanReadFile(this->FileName.c_str()) ) 
   {
       Reader = JPEGReader; 
   }
   else
   {
       return false;
   }
   Reader->SetFileName(this->FileName.c_str());
   Reader->SetFileDimensionality(3);
   return true;
}

void ViewerBase::WriteFile( const std::string& filename ) {
    vtkWUWriter* Writer = vtkWUWriter::New();
    Writer->SetFileName(filename.c_str());
    Writer->SetInput(Reader->GetOutput());
    Writer->SetFileDimensionality(3);
    Writer->Write();
}
