#ifndef __vtkFLTKInstantiator_h
#define __vtkFLTKInstantiator_h

#include "vtkInstantiator.h"
#include "vtkFLTKConfigure.h"

class vtkFLTKInstantiatorInitialize;

class VTK_FLTK_EXPORT vtkFLTKInstantiator
{
  friend class vtkFLTKInstantiatorInitialize;

  static void ClassInitialize();
  static void ClassFinalize();

  static vtkObject* Create_vtkFLTKObjectFactory();
  static vtkObject* Create_vtkFLTKOpenGLRenderWindow();
  static vtkObject* Create_vtkFLTKRenderWindowInteractor();
};

class VTK_FLTK_EXPORT vtkFLTKInstantiatorInitialize
{
public:
  vtkFLTKInstantiatorInitialize();
  ~vtkFLTKInstantiatorInitialize();
private:
  static unsigned int Count;
};

static vtkFLTKInstantiatorInitialize vtkFLTKInstantiatorInitializer;

#endif
