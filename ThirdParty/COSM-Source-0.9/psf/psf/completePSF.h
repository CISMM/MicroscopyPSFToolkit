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
                                                                                
#ifndef _COMPLETE_PSF_H
#define _COMPLETE_PSF_H

#include "psf/radialPSF.h"
#include "psf/psfUser.h"
#include <blitz/array.h>

using namespace blitz;

namespace cosm {

enum PsfType {
    OPTICAL_SECTIONING_WIDEFIELD = 0x0,
    OPTICAL_SECTIONING_2_PHOTON = 0x1,
    CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE = 0x2,
    CONFOCAL_ROTATING_DISK_LINE_APERTURE = 0x3,
	DIC = 0x4,
	DIC_2D = 0x5
};

template<typename T>
class CompletePSF  {

  public:
  
    CompletePSF( 
        RadialPSF<T>* radialPSF,
        PsfType type
    ) : radialPSF_(radialPSF), type_(type) {};

    virtual ~CompletePSF(){};

    virtual void rotate( bool exact = false, EvalType type = MAGNITUDE, PsfUser* user = NULL );
    virtual void rotateAndSum( bool exact = false, EvalType type = MAGNITUDE, PsfUser* user = NULL );
    void sumY( void );
	virtual void rotateXY( double angle ) {};

    Array<T,3> psf() { return psf_; };
	Array<T,3> psfReal() { return psfReal_; };
    Array<T,3> psfImag() { return psfImag_; };

 protected:

    // not allowed
    CompletePSF(CompletePSF<T>&);
    CompletePSF& operator=(CompletePSF<T>&);

  protected:

    RadialPSF<T>* radialPSF_;
    Array<T,3> psf_;
    Array<T,3> psfReal_;
    Array<T,3> psfImag_;
    PsfType type_;

};

};

#include "psf/completePSF.cxx"

#endif // _COMPLETE_PSF_H;
