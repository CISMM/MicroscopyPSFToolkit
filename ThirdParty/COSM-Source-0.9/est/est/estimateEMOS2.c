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

#include "estimateEMOS2.h"
#include "arrayManip.h"
#include "RectDomainIter.h"

namespace cosm {

// class constructor
template<typename T, int N>
EstimateEMOS2<T,N>::EstimateEMOS2(
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
    prev_.resize(extent);
    extent(0) *= 2;
    s2_.resize(extent);
    extent(N-1) = extent(N-1)/2+1;
    estF_.resize(extent);
    sF_.resize(extent);

    // resize the psfs and compute the Fourier transform (OTF)
    for ( int m = 0; m < psfs.length(0); m++ ) 
    {
        psfsF_(m).resize(extent);
        Array<T,N> psfResized(s2_.extent()); 
        padCenter(psfs(m), psfResized); 
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
    // compute the scaling factors
    scale_ = T(1.0)/strata.length(0);
}

// Function to specify the psf for new estimation
template<typename T, int N>
void EstimateEMOS2<T,N>::setPSF(
    Array<T,N>& psf
) {
    // not implemented
}

// Function to specify the image for new estimation
// NOTE: new img has to be the same size as previous img
template<typename T, int N>
void EstimateEMOS2<T,N>::setImage(
    Array<T,N>& img
) {
    EstimateIterative<T,N>::setImage(img);
}

// Function for single iteration of the em algorithm
template<typename T, int N>
void EstimateEMOS2<T,N>::iterate(
){
    int l, m;
    RectDomain<N> rect(this->est_.lbound(), this->est_.ubound());
    this->old_ = this->est_;
    int length = strata_.length(0);
    int size = length - 2;
    for ( l = 0; l < length; l++ ) 
    {
        estF_ = 0;
        prev_ = this->est_;
        // get image prediction which is a convolution of est with 
        // interpolation of psfs (multiplication and in Fourier domain)
        for ( m = 0; m < size; m++ ) 
        {
            s1_ = 0;
            if ( m == 0 )
            {
                s1_(strata_(m)) = this->est_(strata_(m));
            }
            else if ( m == size - 1 )
            {
                s1_(strata_(m+2)) = this->est_(strata_(m+2));
            }
            s1_(strata_(m+1)) = this->est_(strata_(m+1));

            s2_ = s1_;
            multiplyStratum(strata_(m+1), s1_, a_, true);
            multiplyStratum(strata_(m+1), s2_, a_, false);

            mirror(s1_, s_);
            fftw_.plan(s_, sF_);
            fftw_.execute();
            estF_ += sF_ * psfsF_(m);

            mirror(s2_, s_);
            fftw_.plan(s_, sF_);
            fftw_.execute();
            estF_ += sF_ * psfsF_(m+1);
        }
        // convert back to space domain
        fftw_.plan(estF_, s2_);
        fftw_.execute();
        s2_ /= s2_.size();
        this->est_ = s2_(rect);
        // get the ratio of image and prediction
        s_ = where( this->est_ > epsilon_, this->img_/this->est_, this->img_/epsilon_);
        mirror(s_, s2_);
        fftw_.plan(s2_, estF_);
        fftw_.execute();
	
        // convolve the ratio with psf and multiply with old estimate
        this->est_ = 0;
        for ( m = 0; m < size; m++ ) 
        {
            sF_ = estF_ * conj(psfsF_(m));
            // multiply with old estimate
            fftw_.plan(sF_, s2_);
            fftw_.execute();
            s2_ /= s2_.size();
            s_ = s2_(rect);
            s_ *= prev_;
            this->est_ += s_;
        }
        this->est_ *= scale_;
        this->est_ = where( this->est_ > epsilon_, this->est_, 0);
    }
}

}
