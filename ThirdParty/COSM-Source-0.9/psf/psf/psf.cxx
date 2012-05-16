/****************************************************************************
 * Copyright (c) 2007 Einir Valdimarsson and Chrysanthe Preza
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is :wdistributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 ****************************************************************************/
 /****************************************************************************/
 /*  08/11/09 Added parameters to handle dic */
 /*  garyc shear,bias,amplitude Ratio        */
 /****************************************************************************/

#include "psf/gibsonLaniPsfFunctor.h"
#include "psf/haeberlePsfFunctor.h"
#include "xcosm/completeXcosm.h"
#include <ostream>

using namespace cosm;

const float refractiveIndexOfOil = 1.515; 

template<typename T>
void Psf<T>::parameters(
    const PsfParameters<T>& param
) {
   parameters(
      param.Nxy(),
      param.Nz(),
      param.deltaXY(),
      param.deltaZ(),
      param.ts(),
      param.tid(),
      param.tia(),
      param.tgd(),
      param.tga(),
      param.ns(),
      param.nid(),
      param.nia(),
      param.ngd(),
      param.nga(),
      param.tld(),
      param.tla(),
      param.lm(),
      param.na(),
      param.lambda(),
      param.absError(),
      param.fsize(),
      param.distance(),
      param.magY(),
      param.shear(),
      param.bias(),
      param.amplitudeRatio(),
	  param.rotation()
   );
}

template<typename T>
void Psf<T>::parameters(
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
    T fsize,
    T distance,
    T magY,
    T shear,
    T bias,
    T ampRatio,
	T rotation

) {
    nXY_ = Nxy;		
    nZ_ = Nz;		
    deltaXY_ = deltaXY;	
    deltaZ_ = deltaZ;	
    shear_ = shear;
    bias_ = bias;
    ampRatio_ = ampRatio;
	rotation_ = rotation;

    deltaxyNyq_ = (type_ == OPTICAL_SECTIONING_WIDEFIELD) || (type_ == OPTICAL_SECTIONING_2_PHOTON) ?  lambda/(na*4.0) : lambda/(na*8.0);
    deltazNyq_ = refractiveIndexOfOil*lambda/(na*na); 
	fsize_ = fsize;
	distance_ = distance;
	magY_ = magY;
	lm_ = lm;

    if ( model_ == Psf<T>::MODEL_GIBSON_LANI )
    {
	functor_ = new GibsonLaniPsfFunctor<T> (
            ts, tid, tia, tgd, tga,
	    ns, nid, nia, ngd, nga,
	    tld, tla,
	    lm, na, lambda, absError
        );
    }
    else if ( model_ == Psf<T>::MODEL_HAEBERLE )
    {
        functor_ = new HaeberlePsfFunctor<T>(
            ts, tid, tia, tgd, tga,
	    ns, nid, nia, ngd, nga,
	    tld, tla,
	    lm, na, lambda, absError
        );
    }
    else 
    {
	cout <<"Unknown psf model: " << endl;
    }
}

template<typename T>
void Psf<T>::evaluate(
    unsigned short evalType
){
    std::cout <<"Psf<T>::evaluate; enter" << std::endl;
    radialPSF_ = new RadialPSF<T>(
        nXY_, 
        nZ_,
        deltaXY_, 
        deltaZ_, 
        deltaxyNyq_,
        functor_, 
        functor_->isSymmetric()
    );
    if ( (type_ == OPTICAL_SECTIONING_WIDEFIELD) ||
         (type_ == OPTICAL_SECTIONING_2_PHOTON) )
    {
        completePSF_ = new CompletePSF<T>(radialPSF_, type_);
    } 
    else if ( (type_ == CONFOCAL_ROTATING_DISK_CIRCULAR_APERTURE) ||
              (type_ == CONFOCAL_ROTATING_DISK_LINE_APERTURE) || 
              (type_ == DIC || type_ == DIC_2D) )
    {
        completePSF_ = new CompleteXcosm<T>(radialPSF_, type_, fsize_, distance_, magY_, lm_,shear_,bias_,ampRatio_);
    }

    if ( eval_ == Psf<T>::EVAL_EXACT ) 
    {
	    std::cout <<"exact evaluation" << std::endl;
        if ( (evalType & MAGNITUDE) == MAGNITUDE ) completePSF_->rotate(true, MAGNITUDE, user_);
		if ( (evalType & REAL) == REAL ) completePSF_->rotate(true, REAL, user_);
		if ( (evalType & IMAGINARY) == IMAGINARY ) completePSF_->rotate(true, IMAGINARY, user_);
    } 
    else if ( eval_ == Psf<T>::EVAL_INTERPOLATION )
    {
        std::cout <<"interpolation evaluation" << std::endl;
	    radialPSF_->evaluate(user_);
		if ( (evalType & MAGNITUDE) == MAGNITUDE ) completePSF_->rotateAndSum(false, MAGNITUDE, user_);
		if ( (evalType & REAL) == REAL ) completePSF_->rotateAndSum(false, REAL, user_);
		if ( (evalType & IMAGINARY) == IMAGINARY ) completePSF_->rotateAndSum(false, IMAGINARY, user_);
        //completePSF_->rotateAndSum(true, user_);
    }
    else 
    {
	std::cout <<"unknown evaluation method: "<< std::endl;
	return;
    }
	std::cout <<"Psf<T>::evaluate; exit" << std::endl;
}

template<typename T>
void Psf<T>::sumY(
    void
) {
    completePSF_->sumY();
}
