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

#include "estimateEMSV.h"
#include "arrayManip.h"
#include "RectDomainIter.h"

namespace cosm {

// class constructor
template<typename T, int N>
EstimateEMSV<T,N>::EstimateEMSV(
    Array< RectDomain<N>, 1>& strata,
    Array< Array<T,N>, 1>& psfs,
    Array<T,N>& img,
    int iterations,
    T epsilon

) : 
    EstimateIterative<T,N>(psfs(0), img, iterations),
    epsilon_(epsilon), strata_(strata)
{
    a_.resize(this->img_.extent(0));
    strata_.resize(strata.length(0));
    psfsF_.resize(psfs.length(0));
    TinyVector<int,N> extent(this->img_.extent());
    s_.resize(extent);
    s1_.resize(extent);
    extent(N-1) = extent(N-1)/2+1;
    estF_.resize(extent);
    sF_.resize(extent);

    // resize the psfs and compute the Fourier transform (OTF)
    for ( int m = 0; m < psfs.length(0); m++ ) 
    {
		Array<T,N> psfResized(this->img_.extent()); 
        padCenter(psfs(m), psfResized); 
		psfsF_(m).resize(extent);
        psfsF_(m) = forwardFFT(psfResized);
    }

    // compute the interpolation constants
    a_ = 0;
    int index = strata_(0).lbound(int(0));
    for ( int m = 0; m < strata.length(0); m++ ) 
    {
        int l = strata_(m).lbound(int(0));
        int u = strata_(m).ubound(int(0));
	int size = u - l + 1;
	for ( int j = 0; j < size; j++ ) 
	{ 
	    a_(index) = T(1) - T(j)/T(size);
	    index++;
	}
    }
    // std::cout <<"a_: "<< a_ << std::endl;
    // compute the scaling factors
    scale_ = T(1.0)/strata.length(0);
}

// Function to specify the psf for new estimation
template<typename T, int N>
void EstimateEMSV<T,N>::setPSF(
    Array<T,N>& psf
) {
    // not implemented
}

// Function to specify the image for new estimation
// NOTE: new img has to be the same size as previous img
template<typename T, int N>
void EstimateEMSV<T,N>::setImage(
    Array<T,N>& img
) {
    EstimateIterative<T,N>::setImage(img);
}

// Function for single iteration of the em algorithm
template<typename T, int N>
void EstimateEMSV<T,N>::iterate(
){
    int m;
    // get convolution of est with interpolation of psfs 
    // (multiplication and in Fourier domain)
    this->old_ = this->est_;
    estF_ = 0;
    int size = strata_.length(0) - 2;
    for ( m = 0; m < size; m++ ) 
    {
	s_ = 0;
        if ( m == 0 )
        {
	    s_(strata_(m)) = this->est_(strata_(m));
        }
        else if ( m == size - 1 )
        {
	    s_(strata_(m+2)) = this->est_(strata_(m+2));
        }
	s_(strata_(m+1)) = this->est_(strata_(m+1));

        s1_ = s_;
        multiplyStratum(strata_(m+1), s_, a_, true);
        multiplyStratum(strata_(m+1), s1_, a_, false);

        fftw_.plan(s_, sF_);
        fftw_.execute();
	estF_ += sF_ * psfsF_(m);

        fftw_.plan(s1_, sF_);
        fftw_.execute();
	estF_ += sF_ * psfsF_(m+1);
    }
    // convert back to time domain
    fftw_.plan(estF_, this->est_);
    fftw_.execute();
    this->est_ /= this->est_.size();
    // get the ratio of image and convolution
    s_ = where( this->est_ > epsilon_, this->img_/this->est_, this->img_/epsilon_);
    fftw_.plan(s_, estF_);
    fftw_.execute();

    // update estimate by convolving psf with ratio (do it in Fourier domain)
    // and interpolate 
    this->est_ = 0;
    for ( m = 0; m < size; m++ ) 
    {
        sF_ = conj(psfsF_(m)) * estF_; 
        fftw_.plan(sF_, s1_);
        fftw_.execute();
        s1_ /= s1_.size();
	multiplyStratum(strata_(m+1), s1_, a_, true);

        sF_ = conj(psfsF_(m+1)) * estF_; 
        fftw_.plan(sF_, s_);
        fftw_.execute();
        s_ /= s_.size();
	multiplyStratum(strata_(m+1), s_, a_, false);
        s1_(strata_(m+1)) += s_(strata_(m+1));

        s_ = 0;
        if ( m == 0 )
        {
            s_(strata_(m)) = this->old_(strata_(m));
        } 
        else if ( m == size-1 )
        {
            s_(strata_(m+2)) = this->old_(strata_(m+2));
        }
        s_(strata_(m+1)) = this->old_(strata_(m+1));
       
	this->est_ += s_ * s1_;
    }
//    this->est_ *= scale_;
    this->est_ = where( this->est_ > epsilon_, this->est_, 0);
}

}
