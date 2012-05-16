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

#ifndef _ESTIMATE_MAP_H
#define _ESTIMATE_MAP_H

#include "est/estimateLinear.h"
#include <blitz/tinyvec-et.h>

#if WIN32
#define M_PI 3.14159265358979323846
#endif

namespace cosm {

template<typename T, int N>
class EstimateMAP : public EstimateLinear<T,N> {

  public:
    
    // Constructor
    EstimateMAP(
        Array<T,N>& psf,
        Array<T,N>& img,
	T alpha
    ) : EstimateLinear<T,N>(psf, img, alpha) { };

    // Destructor
    virtual ~EstimateMAP() {};

    // Calculate estimate
    virtual int run() 
    {
        TinyVector<double,N> w0 = (2.0*M_PI)/this->psfF_.length();
        ArrayIterator<std::complex<T>,N> iter = this->estF_.begin(), 
				 	   end = this->estF_.end();
        while ( iter != end )
        {
            TinyVector<int,N> pos = iter.position();
            if ( real(this->imgF_(pos)) != 0 || imag(this->imgF_(pos)) != 0 )
            {
                this->estF_(pos) = (this->imgF_(pos)*conj(this->psfF_(pos)))/
	             (this->psfF_(pos)*conj(this->psfF_(pos)) + (T(2) * this->val_) * T(sum(pos * w0)));
            } else {
            	this->estF_(pos) = 0;
	    }
            ++iter;
        }
        this->est_ = inverseFFT(this->estF_);	
        this->est_ = where( this->est_ < 0, 0, this->est_);
	return 0;
    };

    void setAlpha( T alpha ) { this->val_ = alpha; };

    T getAlpha( ) { return this->val_; };

  private:

    // not allowed
    EstimateMAP( EstimateMAP<T,N>& );
    void operator=( EstimateMAP<T,N>& );

};

}

#endif // _ESTIMATE_MAP_H
