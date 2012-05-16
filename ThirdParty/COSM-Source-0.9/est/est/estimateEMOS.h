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
#ifndef _ESTIMATE_EMOS_H
#define _ESTIMATE_EMOS_H

#include "estimateIterative.h"
#include <complex>
#include "fftwInterface.h"

namespace cosm {

template<typename T, int N>
class EstimateEMOS : public EstimateIterative<T,N> {

  public:

    // class constructor
    EstimateEMOS(
	Array< RectDomain<N>, 1>& strata,
	Array< Array<T,N>, 1>& psf,
	Array<T,N>& img,
	int iterations,
	T epsilon = 1E-4
    );

    // class destructor
    virtual ~EstimateEMOS() { };

    // Function to specify the raw image
    virtual void setImage( 
	Array<T,N>& img 
    );

  protected:

    // Function for single iteration of the algorithm
    virtual void iterate();

  private:

    // not allowed
    EstimateEMOS( EstimateEMOS<T,N>& );
    void operator=( EstimateEMOS<T,N>& );

    // Function to specify the psf
    virtual void setPSF( 
	Array<T,N>& psf 
    );
 
  protected:
    
    T epsilon_;				
    fftwInterface<T,N> fftw_;
    Array<Array<std::complex<T>,N>, 1> psfsF_;    
    Array<T,N> prev_;		
    Array<T,N> s_;		
    Array<T,N> s1_;		
    Array<T,N> r_;
    Array<std::complex<T>,N> estF_;  
    Array<std::complex<T>,N> sF_;  
    Array<std::complex<T>,N> rF_;
    Array<RectDomain<N>, 1> strata_;
    Array<RectDomain<N>, 1> subset_;
    Array<int,1> subsetToStrata_;
    Array<T,1> a_;
    T scale_;
};

}

#include "estimateEMOS.c"

#endif  // _ESTIMATE_EMOS_H
