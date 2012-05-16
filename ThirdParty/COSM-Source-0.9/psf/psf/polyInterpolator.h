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

/*
 * Given arrays xa[1..n] and ya[1..n], and given a value x, this routine 
 * returns a value y, and an error estimate dy. If P(x) is the polynomial 
 * of degree N  1 such that P(xai) = yai, i = 1, . . . , n, then the 
 * returned value y = P(x).
 */
#ifndef _POLY_INTERPOLATOR_H
#define _POLY_INTERPOLATOR_H

#include "interpolator.h"
#include <stdlib.h>

namespace cosm {

template<typename T>
class PolyInterpolator : public Interpolator<T> {

  public:

    PolyInterpolator( T x[], T y[], int n ) 
	: Interpolator<T>(x, y, n) 
    {
        c_ = new T[n];
        d_ = new T[n];
    };

    ~PolyInterpolator() 
    {
	delete[] c_;
	delete[] d_;
    };

    // Returns interpolation value
    virtual T operator()( T x )
    {
        int ns = 0;
        T den, dift, ho, hp, w;
	T res;
        T dif = abs(x-this->x_[0]);
        for ( int i = 0; i < this->n_; i++ ) { 
	    // Here we nd the index ns of the closest table entry,
            if ( (dift = abs(x-this->x_[i])) < dif ) 
	    {
                ns = i;
                dif = dift;
            }
	    // and initialize the tableau of cs and ds.
	    c_[i] = this->y_[i]; 
	    d_[i] = this->y_[i];
        }
	// This is the initial approximation.
        res = this->y_[ns--]; 
        // For each column of the tableau, we loop over the current cs and ds 
        // and update them. 
        for ( int m = 1; m < this->n_; m++ ) { 
	    for ( int i = 0; i < this->n_-m; i++ ) { 
	        ho = this->x_[i] - x;
	        hp = this->x_[i+m] - x;
	        w = c_[i+1] - d_[i];
	        // This error can occur only if two input xas are 
	        // (to within roundo) identical.
	        if ( (den = ho-hp) == 0.0 ) 
		{
		    this->ier_ = true;
		    return 0;
		}
	        den = w/den;
		// here cs and ds are updated
	        c_[i] = ho*den;
	        d_[i] = hp*den; 
	    }
            // After each column in the tableau is completed, we decide which 
	    // correction, c or d, we want to add to our accumulating value,
	    //  i.e., which path to take through the tableauforking up or down. 
	    // We do this in such a way as to take the most straight line route 
	    // through the tableau to its apex, updating ns accordingly to keep 
	    // track of where we are. This route keeps the partial 
	    // approximations centered (insofar as possible) on the target x. 
	    // The last dy added is thus the error indication.
	    this->error_ = 2*ns < (this->n_-m) ? c_[ns+1] : d_[ns--];
	    res += this->error_; 
        }
	this->ier_ = false;
	return res;
    };

    virtual void setValues( T x[], T y[], int n) 
    { 
	if ( n > this->n_ ) 
	{
	    delete[] c_;
	    delete[] d_;
	    c_ = new T[n];
	    d_ = new T[n];
	}
	Interpolator<T>::setValues(x,y,n);
    };

  protected:

    // not allowed
    PolyInterpolator(PolyInterpolator<T>&);
    PolyInterpolator& operator=(PolyInterpolator<T>&);

  protected:

    T* c_;
    T* d_;

};

}

#endif // _INTERPOLATOR_H
