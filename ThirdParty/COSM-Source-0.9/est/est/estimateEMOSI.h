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
#ifndef _ESTIMATE_EMOSI_H
#define _ESTIMATE_EMOSI_H

#include "estimateIterative.h"
#include <complex>
#include "fftwInterface.h"

namespace cosm {

template<typename T, int N>
class EstimateEMOSI : public EstimateIterative<T,N> {

  public:

    // class constructor
    EstimateEMOSI(
	Array<T,N>& psf,
	Array<T,N>& img,
	int iterations,
	EstimateUser<T,N>* user = NULL,
	int update = -1,
	T epsilon = 1E-4
    );

    // class destructor
    virtual ~EstimateEMOSI() { };

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

    void subset( 
	int l,
	Array<T,2>& img,
	Array<T,2>& est,
	Array<T,2>& s,
	Array<T,2>& old,
	Array<std::complex<T>,2>& psfF
    );

    void subset( 
	int l,
	Array<T,3>& img,
	Array<T,3>& est,
	Array<T,3>& s,
	Array<T,3>& old,
	Array<std::complex<T>,3>& psfF
   );

  private:

    // not allowed
    EstimateEMOSI( EstimateEMOSI<T,N>& );
    void operator=( EstimateEMOSI<T,N>& );
 
  protected:
    
    T epsilon_;				
    fftwInterface<T,N> fftw_;
    fftwInterface<T,N-1> subfftw_;
    Array<std::complex<T>,N> psfF_;     // OTF - Fourier transform of psf
    Array<std::complex<T>,N> sF_;      
    Array<T,N> s_;      
    Array<std::complex<T>,N-1> rF_;      
};

}

#include "estimateEMOSI.c"

#endif  // _ESTIMATE_EMOSI_H
