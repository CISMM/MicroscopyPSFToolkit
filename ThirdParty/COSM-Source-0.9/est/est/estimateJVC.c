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

#include "est/estimateJVC.h"

namespace cosm {

// class constructor
template<typename T, int N>
EstimateJVC<T,N>::EstimateJVC(
    Array<T,N>& psf,
    Array<T,N>& img,
    int iterations,
    T epsilon
) : 
    EstimateIterative<T,N>(psf, img, iterations),
    A_(0),
    epsilon_(epsilon)
{
    TinyVector<int,N> extent(this->img_.extent());
    g_.resize(extent);
    extent(N-1) = extent(N-1)/2+1;
    estF_.resize(extent);
    setPSF(psf);
    A_ = (max)(this->img_)/T(2.0);
    cout <<"A_: "<< A_ << endl;
}

// Function to specify the psf for new estimation
template<typename T, int N>
void EstimateJVC<T,N>::setPSF(
    Array<T,N>& psf
) {
    EstimateIterative<T,N>::setPSF(psf);
    TinyVector<int,N> extent(this->img_.extent());
    extent(N-1) = extent(N-1)/2+1;
    psfF_.resize(extent);

    Array<T,N> psfResized(this->img_.extent()); 
    padCenter(this->psf_, psfResized); 
    psfF_ = forwardFFT(psfResized);
    //psfF_ /= sum(psfF_);

    this->est_ = this->img_;  		// default initial estimation guess	
}

// Function to specify the image for new estimation
template<typename T, int N>
void EstimateJVC<T,N>::setImage(
    Array<T,N>& img
) {
    bool resizeFlag =  !equal(img.shape(),this->img_.shape()) ? true : false;
    EstimateIterative<T,N>::setImage(img);
    if ( resizeFlag ) {
        TinyVector<int,N> extent(this->img_.extent());
        g_.resize(extent);
        Array<T,N> psfResized(extent); 
        extent(N-1) = extent(N-1)/2+1;
        psfF_.resize(extent);
        estF_.resize(extent);
        padCenter(this->psf_, psfResized); 
        psfF_ = forwardFFT(psfResized);
	//psfF_ /= sum(psfF_);
    }
    A_ = (max)(this->img_)/T(2.0);
    this->est_ = this->img_;       	// default initial estimation guess	
}

// Function for single iteration of the em algorithm
template<typename T, int N>
void EstimateJVC<T,N>::iterate(
){
    // get convolution of est with psf (multiplication in Fourier domain)
    this->old_ = this->est_;
    fftw_.plan(this->est_, estF_);
    fftw_.execute();
    estF_ = estF_ * psfF_;
    // convert back to time domain
    fftw_.plan(estF_, g_);
    fftw_.execute();
    g_ /= g_.size();
    // update the estimate
//    this->est_ = g_;
//    this->est_ = T(1) -  pow2( (this->old_ - A_)/A_);
//    this->est_ = this->img_ - g_;
    this->est_ = this->old_ + (T(1) - pow2( (this->old_ - A_)/A_) ) * ( this->img_ - g_ );
    // apply constraints
    T minVal = T(0);
    T maxVal = T(2)*A_;
    this->est_ = where( this->est_ < minVal, minVal, this->est_);
    this->est_ = where( this->est_ > maxVal, maxVal, this->est_);
    if ( (max)( abs(this->img_ - g_)) < epsilon_ ) 
    {
	abort();
    }
}

}
