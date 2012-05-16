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

#ifndef _ESTIMATE_LINEAR_H
#define _ESTIMATE_LINEAR_H

#include "est/estimate.h"
#include <complex>
#include "blitz/fftwInterface.h"
#include "blitz/RectDomainIter.h"
#include "blitz/arrayManip.h"

namespace cosm {

template<typename T, int N>
class EstimateLinear : public Estimate<T,N> {

  public:
    
    // Class constructor
    EstimateLinear(
        Array<T,N>& psf,
        Array<T,N>& img,
	T val
    ) : Estimate<T,N>(psf, img), val_(val) 
    {
	TinyVector<int, N> extent(this->img_.extent());
	// real fft has different last dimension 
        extent(N-1) = extent(N-1)/2+1;
	// adjust to the correct sizes
	imgF_.resize(extent);
	estF_.resize(extent);
        setPSF(psf);
	imgF_ = forwardFFT(this->img_);
    };

    // Class destructor
    virtual ~EstimateLinear() {};

    // Calculate estimate
    virtual int run() = 0;

    // Function to specify the psf
    virtual void setPSF(
        Array<T,N>& psf
    ) { 
	Estimate<T,N>::setPSF(psf); 
	TinyVector<int, N> extent(this->img_.extent());
	// real fft has different last dimension 
        extent(N-1) = extent(N-1)/2+1;
	// adjust to the correct sizes
	psfF_.resize(extent);
	Array<T,N> psfResized(this->img_.extent()); 
	padCenter(psf, psfResized);
	// take the fft of psf
	psfF_ = forwardFFT(psfResized);
	this->est_ = 0;
        wuDataWrite(psfResized, "psf_lls.wu");
    };

    // Function to specify the raw image
    virtual void setImage(
        Array<T,N>& img
    ) { 
	bool resizeFlag = !equal(img.shape(),this->img_.shape()) ? true : false;
	Estimate<T,N>::setImage(img); 
	if ( resizeFlag ) {
	    TinyVector<int, N> extent(this->img_.extent());
 	    psfF_.resize(imgF_.extent());
	    Array<T,N> psfResized(this->img_.extent()); 
	    padCenter(this->psf_,psfResized);
	    psfF_ = forwardFFT(psfResized);
	    estF_.resize(imgF_.extent());
	    imgF_.resize(extent);
	}
	imgF_ = forwardFFT(this->img_);
	this->est_ = 0;
    };

    // Function to set algorithm value
    void setVal( T val ) { val_ = val; };

    // Function to get algorithm value
    T getVal() { return val_; };

  protected:
    
    Array<std::complex<T>,N> psfF_;     // OTF - Fourier transform of psf
    Array<std::complex<T>,N> imgF_;     // Fourier transform of image
    Array<std::complex<T>,N> estF_;     // Fourier transform of estimate
    T val_;

};

}

#endif // _ESTIMATE_LINEAR_H
