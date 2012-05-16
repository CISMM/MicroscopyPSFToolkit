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

#ifndef _ESTIMATE_EMOS2_BIG_H
#define _ESTIMATE_EMOS2_BIG_H

#include "estimateIterative.h"
#include "estimateIO.h"
#include <complex>
#include "fftwInterface.h"

namespace cosm {

template<typename T, int N>
class EstimateEMOS2big : public EstimateIterative<T,N> {

  public:

    // class constructor
    EstimateEMOS2big(
	Array< RectDomain<N>, 1>& strata,
        const std::string& psfName,
        const std::string& otfName,
        const std::string& suffix,
	Array<T,N>& img,
	int iterations,
        EstimateIO<T,N>* io,
	T epsilon = 1E-4
    );

    // class destructor
    virtual ~EstimateEMOS2big() { };

    // Function to specify the raw image
    virtual void setImage( 
	Array<T,N>& img 
    );

  protected:

    // Function for single iteration of the algorithm
    virtual void iterate();

  private:

    // not allowed
    EstimateEMOS2big( EstimateEMOS2big<T,N>& );
    void operator=( EstimateEMOS2big<T,N>& );

    // Function to specify the psf
    virtual void setPSF( 
	Array<T,N>& psf 
    );
 
  protected:
    
    T epsilon_;				
    Array<RectDomain<N>, 1> strata_;
    std::vector<std::string> otfName_;
    fftwInterface<T,N> fftw_;
    Array<T,N> prev_;		
    Array<T,N> s_;		
    Array<T,N> s2_;		
    Array<std::complex<T>,N> psfF_;  
    Array<std::complex<T>,N> estF_;  
    Array<std::complex<T>,N> sF_;  
    Array<T,1> a_;
    T scale_;
    EstimateIO<T,N>* io_;
};

}

#include "estimateEMOS2big.c"

#endif  // _ESTIMATE_EMOS2_BIG_H
