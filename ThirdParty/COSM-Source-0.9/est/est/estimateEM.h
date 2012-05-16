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
#ifndef _ESTIMATE_EM_H
#define _ESTIMATE_EM_H

#include "est/estimateIterative.h"
#include <complex>
#include "blitz/fftwInterface.h"

namespace cosm {

template<typename T, int N>
class EstimateEM : public EstimateIterative<T,N> {

  public:

    // class constructor
    EstimateEM(
	Array<T,N>& psf,
	Array<T,N>& img,
	int iterations,
        EstimatePenalty<T,N>* penalty = NULL,
	T epsilon = 1E-4
    );

    // class destructor
    virtual ~EstimateEM() { };

    // Function to specify the psf
    virtual void setPSF( 
	Array<T,N>& psf 
    );

    // Function to specify the raw image
    virtual void setImage( 
	Array<T,N>& img 
    );

  protected:

    // Function for single iteration of the em algorithm
    virtual void iterate();

  private:

    // not allowed
    EstimateEM( EstimateEM<T,N>& );
    void operator=( EstimateEM<T,N>& );
 
  protected:
    
    T epsilon_;				
    fftwInterface<T,N> fftw_;
    Array<std::complex<T>,N> psfF_;     // OTF - Fourier transform of psf
    Array<std::complex<T>,N> estF_;     // Fourier transform of estimate

};

}

#include "estimateEM.c"

#endif  // _ESTIMATE_EM_H
