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
                                                                                
#ifndef _COMPLETE_XCOSM_H
#define _COMPLETE_XCOSM_H

#include "psf/completePSF.h"
#include "psf/radialPSF.h"
#include "psf/psfUser.h"
#include <blitz/array.h>

namespace cosm {

template<typename T>
class CompleteXcosm : public CompletePSF<T> {

  public:
  
    CompleteXcosm( 
        RadialPSF<T>* radialPSF,
        PsfType type,
        T fsize = 0,
        T distance = 0,
        T magY = 1, 
        T lm_ = 1,
        T shear=0,
        T bias=0,
        T ampRatio=0
    ) : CompletePSF<T>(radialPSF, type), fsize_(fsize), distance_(distance), magY_(magY),shear_(shear),bias_(bias),ampRatio_(ampRatio) {};

    ~CompleteXcosm(){};

    virtual void rotate( bool exact = false, PsfUser* user = NULL );
    virtual void rotateAndSum( bool exact = false, PsfUser* user = NULL );
	virtual void rotateXY( double angle );

 protected:
 
    void rotateDIC( PsfUser* user );

    // not allowed
    CompleteXcosm(CompleteXcosm<T>&);
    CompleteXcosm& operator=(CompleteXcosm<T>&);

protected:
    T fsize_;
    T distance_;
	T magY_;
	T lm_;
    T shear_;
    T bias_;
    T ampRatio_;
    blitz::Array<float,2> radial_;
	blitz::Array<float,2> radialRe_;
	blitz::Array<float,2> radialIm_;
    blitz::Array<float,3> complete_;
	blitz::Array<float,3> completeRe_;
    blitz::Array<float,3> completeIm_;


};

};

#include "xcosm/completeXcosm.cxx"

#endif // _COMPLETE_XCOSM_H;
