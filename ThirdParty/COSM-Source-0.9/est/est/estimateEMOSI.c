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

#include "estimateEMOSI.h"
#include "arrayManip.h"
#include "RectDomainIter.h"

namespace cosm {

// class constructor
template<typename T, int N>
EstimateEMOSI<T,N>::EstimateEMOSI(
    Array<T,N>& psf,
    Array<T,N>& img,
    int iterations,
    EstimateUser<T,N>* user,
    int update,
    T epsilon
) : 
    EstimateIterative<T,N>(psf, img, iterations, user, update ),
    epsilon_(epsilon)
{
    TinyVector<int,N-1> subextent;
    for ( int i = 0; i < N-1; i++ )
    {
	subextent(i) = this->img_.length(i+1);
    }
    subextent(N-2) = subextent(N-2)/2+1;
    rF_.resize(subextent);

    TinyVector<int,N> extent(this->img_.shape());
    s_.resize(extent);
    extent(N-1) = extent(N-1)/2+1;
    psfF_.resize(extent);
    sF_.resize(extent);
    if ( resizeFlag ) {
        Array<T,N> psfResized(this->img_.shape()); 
        padCenter(this->psf_, psfResized); 
        psfF_ = forwardFFT(psfResized);
    } else {
        psfF_ = forwardFFT(psf);
    }
    if ( this->user_ != NULL ) 
    {
	this->user_->initialize(this->est_);     // user supplied initial guess
    } 
    else 
    {
        this->est_ = 1;  		// default initial estimation guess	
    }
}

// Function to specify the psf for new estimation
template<typename T, int N>
void EstimateEMOSI<T,N>::setPSF(
    Array<T,N>& psf
) {
    bool resizeFlag = !equal(this->psf_.shape(), this->img_.shape());
    EstimateIterative<T,N>::setPSF(psf);
    if ( resizeFlag ) {
        TinyVector<int,N> extent(this->img_.shape());
        extent(N-1) = extent(N-1)/2+1;
        Array<T,N> psfResized(this->img_.shape()); 
        padCenter(this->psf_, psfResized); 
        psfF_.resize(extent);
        psfF_ = forwardFFT(psfResized);
	psfF_ /= sum(psfF_);
    } else {
        psfF_ = forwardFFT(psf);
    }
    if ( this->user_ != NULL ) 
    {
	this->user_->initialize(this->est_);     // user supplied initial guess
    } 
    else 
    {
        this->est_ = 1;  		// default initial estimation guess	
    }
}

// Function to specify the image for new estimation
template<typename T, int N>
void EstimateEMOSI<T,N>::setImage(
    Array<T,N>& img
) {
    bool resizeFlag =  !equal(img.shape(),this->img_.shape()) ? true : false;
    EstimateIterative<T,N>::setImage(img);
    if ( resizeFlag ) {
        TinyVector<int,N-1> subextent;
        for ( int i = 0; i < N-1; i++ )
        {
	    subextent(i) = this->img_.length(i+1);
        }
        subextent(N-2) = subextent(N-2)/2+1;
        rF_.resize(subextent);

        TinyVector<int,N> extent(this->img_.shape());
        s_.resize(extent);
        Array<T,N> psfResized(extent); 
        extent(N-1) = extent(N-1)/2+1;
        psfF_.resize(extent);
        sF_.resize(extent);
        padCenter(this->psf_, psfResized); 
        psfF_ = forwardFFT(psfResized);
	psfF_ /= sum(psfF_);
    }
    if ( this->user_ != NULL ) 
    {
	this->user_->initialize(this->est_);     // user supplied initial guess
    } 
    else 
    {
        this->est_ = 1;  		// default initial estimation guess	
    }
}

// Function for single iteration of the em algorithm
template<typename T, int N>
void EstimateEMOSI<T,N>::iterate(
){
    // get convolution of est with psf (multiplication in Fourier domain)
    this->old_ = this->est_;
    int L = this->est_.length(0);
   
    for ( int l = 0; l < L; l++ ) 
    {
        fftw_.plan(this->est_, sF_);
        fftw_.execute();
        sF_ = sF_ * psfF_;
        // convert back to time domain
        fftw_.plan(sF_, s_);
        fftw_.execute();
        s_ /= s_.size();
	subset(l, this->img_, this->est_, s_, this->old_, psfF_);
	cout <<"this->est_: "<< this->est_ <<endl;
    }
}

template<typename T, int N>
void EstimateEMOSI<T,N>::subset(
    int l,
    Array<T,2>& img,
    Array<T,2>& est,
    Array<T,2>& s,
    Array<T,2>& old,
    Array<std::complex<T>,2>& psfF
){
    Range all = Range::all();
    Array<T, N-1> p = this->psf_(l, all);
    Array<T,N-1> r = s(l,all);
    cout <<"s: "<< r << endl;
    Array<T,N-1> i = img(l,all); 
    cout <<"i: "<< i << endl;
    // get the ratio of image and convolution
    r = where( r > epsilon_, i/r, i/epsilon_);
    cout <<"ratio: "<< r << endl;
    cout <<"p: "<< p << endl;
    // update estimate by convolving psf with ratio 
    Array<std::complex<T>, N-1> pF(rF_.shape());
    subfftw_.plan(p, pF);
    subfftw_.execute();
    subfftw_.plan(r, rF_);
    subfftw_.execute();
    rF_ = conj(pF) * rF_; 
    subfftw_.plan(rF_, r);
    subfftw_.execute();
    r /= r.size();
    cout <<"r: "<< r << endl;
    // multiply with old estimate
    r *= old(l,all);
    r = where( r > epsilon_, r, 0);
    est(l,all) = r;
}

template<typename T, int N>
void EstimateEMOSI<T,N>::subset(
    int l,
    Array<T,3>& img,
    Array<T,3>& est,
    Array<T,3>& s,
    Array<T,3>& old,
    Array<std::complex<T>,3>& psfF
){
    Range all = Range::all();
    Array<T,N-1> r = s(l,all,all);
    Array<std::complex<T>, N-1> pF = psfF(l,all,all);
    Array<T,N-1> i = img(l,all,all); 
    // get the ratio of image and convolution
    r = where( r > epsilon_, i/r, i/epsilon_);
    // update estimate by convolving psf with ratio 
    subfftw_.plan(r, rF_);
    subfftw_.execute();
    rF_ = conj(pF) * rF_; 
    subfftw_.plan(rF_, r);
    subfftw_.execute();
    r /= r.size();
    // multiply with old estimate
    r *= old(l,all,all);
    r = where( r > epsilon_, r, 0);
    est(l,all,all) = r;
}

}
