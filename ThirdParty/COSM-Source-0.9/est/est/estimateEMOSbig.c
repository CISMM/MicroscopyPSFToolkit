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

#include "estimateEMOSbig.h"
#include "arrayManip.h"
#include "RectDomainIter.h"

namespace cosm {

// class constructor
template<typename T, int N>
EstimateEMOSbig<T,N>::EstimateEMOSbig(
    Array< RectDomain<N>, 1>& strata,
    const std::string& psfName,
    const std::string& otfName,
    const std::string& suffix,
    Array<T,N>& img,
    int iterations,
    EstimateIO<T,N>* io,
    T epsilon
) : 
    EstimateIterative<T,N>(img, iterations),
    epsilon_(epsilon), 
    strata_(strata), 
    otfName_(strata.length(0)+1),
    io_(io)
{
    a_.resize(this->img_.extent(0));
    strata_.resize(strata.length(0));
    TinyVector<int,N> extent(this->img_.extent());
    s_.resize(extent);
    prev_.resize(extent);
    extent(N-1) = extent(N-1)/2+1;
    psfF_.resize(extent);
    estF_.resize(extent);
    sF_.resize(extent);

    // resize the psfs and compute the Fourier transform (OTF)
    for ( int m = 0; m <= strata.length(0); m++ ) 
    {

        std::ostringstream psfm;
        psfm << psfName << m << "." << suffix;
        io_->ReadData(this->psf_, psfm.str());
        Array<T,N> psfResized(this->img_.extent()); 
        padCenter(this->psf_, psfResized); 
        psfF_ = forwardFFT(psfResized);
        std::ostringstream otfm;
        otfm << otfName << m << "." << suffix;
        otfName_[m] = otfm.str();
        io_->WriteData(psfF_, otfm.str());
    }

    // compute the interpolation constants
    int index = 0;
    for ( int m = 0; m < strata.length(0); m++ ) 
    {
        int l = strata_(m).lbound(int(0));
        int u = strata_(m).ubound(int(0));
	int size = u - l ;
	for ( int j = 0; j <= size; j++ ) 
	{ 
	    a_(index) = T(1) - T(j)/T(size);
	    index++;
	}
	
    }

    // compute the scaling factors
    scale_ = T(1.0)/strata.length(0);

    // estimate initialization
}

// Function to specify the psf for new estimation
template<typename T, int N>
void EstimateEMOSbig<T,N>::setPSF(
    Array<T,N>& psf
) {
    // not implemented
}

// Function to specify the image for new estimation
// NOTE: new img has to be the same size as previous img
template<typename T, int N>
void EstimateEMOSbig<T,N>::setImage(
    Array<T,N>& img
) {
    EstimateIterative<T,N>::setImage(img);
}

// Function for single iteration of the em algorithm
template<typename T, int N>
void EstimateEMOSbig<T,N>::iterate(
){
    int l, m;
    this->old_ = this->est_;
    int size = strata_.length(0);
    for ( l = 0; l < size; l++ ) 
    {
        estF_ = 0;
        prev_ = this->est_;
	io_->ReadData(psfF_, otfName_[0]);
        // get image prediction which is a convolution of est with 
	// interpolation of psfs (multiplication and in Fourier domain)
        for ( m = 0; m < size; m++ ) 
        {
	    s_ = 0;
	    s_(strata_(m)) = this->est_(strata_(m));
            multiplyStratum(strata_(m), s_, a_, true);
            fftw_.plan(s_, sF_);
            fftw_.execute();
	    estF_ += sF_ * psfF_;
            s_ = 0;
            s_(strata_(m)) = this->est_(strata_(m));
            multiplyStratum(strata_(m), s_, a_, false);
            fftw_.plan(s_, sF_);
            fftw_.execute();
            io_->ReadData(psfF_, otfName_[m+1]);
	    estF_ += sF_ * psfF_;
	}
        // convert back to space domain
        fftw_.plan(estF_, this->est_);
        fftw_.execute();
        this->est_ /= this->est_.size();

        // get the ratio of image and prediction
        s_ = where( this->est_ > epsilon_, this->img_/this->est_, this->img_/epsilon_);
        fftw_.plan(s_, estF_);
        fftw_.execute();
	
	// convolve the ratio with psf and multiply with old estimate
	this->est_ = 0;
        for ( m = 0; m < size; m++ ) 
	{
            io_->ReadData(psfF_, otfName_[m]);
	    sF_ = estF_ * conj(psfF_);
            // multiply with old estimate
	    fftw_.plan(sF_, s_);
	    fftw_.execute();
	    s_ /= s_.size();
            s_ *= prev_;
	    this->est_ += s_;
	}
        this->est_ *= scale_;
        this->est_ = where( this->est_ > epsilon_, this->est_, 0);
    }
}

}
