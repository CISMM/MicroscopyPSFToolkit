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

// Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine 
// returns a value of y and an accuracy estimate. The value returned is 
// that of the diagonal rational function, evaluated at x, which passes 
// through the n points (xai, yai), i = 1...n.

#ifndef _RAT_INTERPOLATOR_H
#define _RAT_INTERPOLATOR_H

#include "interpolator.h"

namespace cosm {

//A small number.
#define TINY 1.0e-25 

template<typename T>
class RatInterpolator : public Interpolator<T> {

  public:

    RatInterpolator( T x[], T y[], int n ) 
	: Interpolator<T>(x, y, n) 
    {
	c_ = new T[n];
	d_ = new T[n];
    };

    ~RatInterpolator() 
    {
	delete[] c_;
	delete[] d_;
    };

    // Returns interpolation value
    virtual T operator()( T x )
    {
	ier_ = false;
	error_ = 0.0;
        int ns = 0;
        T w,t,h,dd;
        T hh = std::abs(x-x_[0]);
        for ( int i = 0; i < n_; i++ ) 
	{
            h = std::abs(x-x_[i]);
            if ( h == 0.0 ) 
	    {
                return y_[i];
            } 
	    else if ( h < hh ) 
	    {
                ns = i;
                hh = h;
            }
            c_[i] = y_[i];
	    // The TINY part is needed to prevent a rare zero-over-zero 
	    // condition. 
            d_[i] = y_[i] + TINY; 
        }
        T y = y_[ns--];
        for ( int m = 1; m <= n_-1; m++ ) 
	{
	    for ( int i = 0; i < n_-m; i++ ) 
	    {
	        w = c_[i+1]-d_[i];
	        // h will never be zero, since this was tested in the 
	        // initializing loop
                h = x_[i+m] - x; 
	        t = (x_[i]-x)*d_[i]/h;
	        dd = t - c_[i+1];
	        if ( dd == 0.0 ) 
		{
	            // This error condition indicates that the interpolating 
		    // function  has a pole at the requested value of x.
		    ier_ = true;
		    return 0;
		}
	        dd = w/dd;
	        d_[i] = c_[i+1] * dd;
	        c_[i] = t*dd;
            }
	    y += (error_ = (2*ns < (n_-m) ? c_[ns+1] : d_[ns--]));
        }
	return y;
    };

    virtual void setValues( T x[], T y[], int n) 
    { 
	if ( n > n_ ) 
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
    RatInterpolator(RatInterpolator<T>&);
    RatInterpolator& operator=(RatInterpolator<T>&);

  protected:

    T* c_;
    T* d_;

};

}

#endif // _RAT_INTERPOLATOR_H
