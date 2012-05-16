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

#ifndef _PSF_H
#define _PSF_H

#include "psf/psfFunctor.h"
#include "psf/radialPSF.h"
#include "psf/completePSF.h"
#include "psf/psfUser.h"
#include "psf/psfParameters.h"
#include <blitz/array.h>

namespace cosm {

template<typename T>
class Psf  {

  public:

    enum Model {
        MODEL_GIBSON_LANI = 0x0,
        MODEL_HAEBERLE = 0x1,
    };
    enum Eval {
        EVAL_INTERPOLATION = 0x0,
        EVAL_EXACT = 0x1
    };

  public:
  
    Psf( Model model = MODEL_GIBSON_LANI,
	 PsfType type = OPTICAL_SECTIONING_WIDEFIELD,
	 Eval eval = EVAL_INTERPOLATION,
	 PsfUser* user = NULL
    ) : nXY_(0), nZ_(0), deltaXY_(0), deltaZ_(0), deltaxyNyq_(0), deltazNyq_(0), shear_(0),bias_(0),ampRatio_(0), rotation_(0), model_(model), type_(type), eval_(eval), user_(user) {};

    ~Psf(){};

    void parameters(
        int Nxy,        // size of the image in x and y
        int Nz,         // size of the image in z
        T deltaXY,      // pixel size in z and y of the image (mm)
        T deltaZ,       // pixel size in z of the image (mm)
        T ts,           // specimen thickness
        T tid,          // immersion thickness design
        T tia,          // immersion thickness actual
        T tgd,          // coverglass thickness design
        T tga,          // coverglass thickness actual
        T ns,           // specimen refractive index
        T nid,          // immersion refractive index design
        T nia,          // immersion refractive index actual
        T ngd,          // coverglass refractive index design
        T nga,          // coverglass refractive index actual
        T tld,          // tube length design
        T tla,          // tube length actual
        T lm,           // lateral magnification
        T na,           // numerical aperture (NA)
        T lambda,
        T absError,
        T fsize = 0,
        T distance = 0,
        T magY = 0,
        T shear=0,
        T bias=0,
        T ampRatio=0,
		T rotation=0
    );

    void parameters(
        const PsfParameters<T>& psfParameters
    );

    void evaluate(
        unsigned short evalType = MAGNITUDE
    );

    void sumY( void );
	void rotateXY( double angle ) { return completePSF_->rotateXY(angle); };

    blitz::Array<T,3> psf() { return completePSF_->psf(); };
    blitz::Array<T,3> psfReal() { return completePSF_->psfReal(); };
	blitz::Array<T,3> psfImag() { return completePSF_->psfImag(); };

 protected:

    // not allowed
    Psf( const Psf&);
    const Psf& operator=( const Psf&);

  protected:

    int nXY_;            // size of the image in x and y
    int nZ_;             // size of the image in z
    T deltaXY_;        // pixel size in z and y of the image (mm)
    T deltaZ_;         // pixel size in z of the image (mm)
    T deltaxyNyq_;     // Nyquist pixel size in x and y
    T deltazNyq_;      // Nyquist pixel size in z
    T fsize_;
    T distance_;
    T magY_;
    T lm_;
    T shear_;
    T bias_;
    T ampRatio_;
	T rotation_;

    Model model_;
    PsfType type_;
    Eval eval_;
    PsfUser* user_;
    RadialPSF<T>* radialPSF_;
    CompletePSF<T>* completePSF_;
    PsfFunctor<T>* functor_;

};

}

#include "psf/psf.cxx"

#endif // _PSF_H;
