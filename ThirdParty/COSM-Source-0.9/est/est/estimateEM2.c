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

#include "est/estimateEM2.h"
#include "blitz/arrayManip.h"
#include "blitz/RectDomainIter.h"
#include <blitz/tinyvec-et.h>

namespace cosm {

// class constructor
template<typename T, int N>
EstimateEM2<T,N>::EstimateEM2(
    Array<T,N>& psf,
    Array<T,N>& img,
    int iterations,
    EstimatePenalty<T,N>* penalty,
    T epsilon
) : 
    EstimateIterative<T,N>(psf, img, iterations, penalty ),
    epsilon_(epsilon)
{
    TinyVector<int,N> extent(this->img_.extent());
    extent(0) = extent(0)*2;
    est2_.resize(extent);
    extent(N-1) = extent(N-1)/2+1;
    estF_.resize(extent);
    setPSF(psf);
}

// Function to specify the psf for new estimation
template<typename T, int N>
void EstimateEM2<T,N>::setPSF(
    Array<T,N>& psf
) {
    EstimateIterative<T,N>::setPSF(psf);
    TinyVector<int,N> extent(this->img_.extent());
    extent(0) = extent(0)*2;
    Array<T,N> psfResized(extent);
    padCenter(this->psf_,psfResized);
    extent(N-1) = extent(N-1)/2+1;
    psfF_.resize(extent);
    psfF_ = forwardFFT(psfResized);
    //psfF_ /= sum(psfF_);
}

// Function to specify the image for new estimation
template<typename T, int N>
void EstimateEM2<T,N>::setImage(
    Array<T,N>& img
) {
    bool resize = !equal(img.shape(),this->img_.shape()) ? true : false;
    EstimateIterative<T,N>::setImage(img);
    if ( resize )
    {
        TinyVector<int,N> extent(this->img_.extent());
        extent(0) = extent(0)*2;
        est2_.resize(extent);
        Array<T,N> psfResized(extent);
        padCenter(this->psf_,psfResized);
        extent(N-1) = extent(N-1)/2+1;
        psfF_.resize(extent);
        estF_.resize(extent);
        psfF_ = forwardFFT(psfResized);
    }
}

// Function for single iteration of the em algorithm
template<typename T, int N>
void EstimateEM2<T,N>::iterate(
) {
    // get convolution of est with psf (multiplication in Fourier domain)
    this->old_ = this->est_;
    mirror(this->est_, est2_);
    fftw_.plan(est2_, estF_);
    fftw_.execute();
    estF_ = estF_ * psfF_;
    // convert back to time domain
    fftw_.plan(estF_, est2_);
    fftw_.execute();
    est2_ /= est2_.size();
    // get the ratio of image and convolution
    RectDomain<N> rect(this->est_.lbound(), this->est_.ubound());
    this->est_ = est2_(rect); 
    this->est_ = where( this->est_ > epsilon_, this->img_/this->est_, this->img_/epsilon_ );
    mirror(this->est_, est2_);
    // update estimate by convolving psf with ratio (do it in Fourier domain)
    fftw_.plan(est2_, estF_);
    fftw_.execute();
    estF_ = conj(psfF_) * estF_; 
    fftw_.plan(estF_, est2_);
    fftw_.execute();
    est2_ /= est2_.size();
    this->est_ = est2_(rect); 
    this->est_ *= this->old_;
    this->est_ = where( this->est_ > epsilon_, this->est_, 0 );
    if ( this->penalty_ != NULL ) 
    {
	this->penalty_->operator()(this->est_);
    }
}

}
