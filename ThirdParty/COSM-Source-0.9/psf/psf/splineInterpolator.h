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

// Given the arrays x[1..n] and y[1..n], which tabulate a function 
// (with the xais in order), and given the array y2_[1..n], which is 
// the output from spline, and given a value of x, this routine 
// returns a cubic-spline interpolated value y.

#ifndef _SPLINE_INTERPOLATOR_H
#define _SPLINE_INTERPOLATOR_H

#include "interpolator.h"

namespace cosm {

template<typename T>
class SplineInterpolator : public Interpolator<T> {

  public:

    SplineInterpolator( T x[], T y[], int n ) 
	: Interpolator<T>(x, y, n) 
    {
	y2_ = new T[n+1];
	u_ = new T[n+1];
	spline();
    };

    ~SplineInterpolator() 
    {
	delete[] y2_;
	delete[] u_;
    };

    // Returns interpolation value
    virtual T operator()( T x )
    {
	int klo,khi,k;
	T h,b,a;
	// We will find the right place in the table by means of
	// bisection. This is optimal if sequential calls to this
	// routine are at random values of x. If sequential calls
	// are in order, and closely spaced, one would do better
	// to store previous values of klo and khi and test if
	// they remain appropriate on the next call.
	klo = 0; 
	khi = n_-1;
	while ( khi-klo > 1 ) {
	    k = (khi+klo) >> 1;
	    if ( x_[k] > x ) 
	    {
		khi = k;
	    }
	    else 
	    {
		klo = k;
	    }
	} 
	// klo and khi now bracket the input value of x.
	h = x_[khi]-x_[klo];
	// The xas must be distinct.
	if ( h == 0.0 ) 
	{
	   ier_ = true;
	   return 0;
	}
	a = (x_[khi]-x)/h;
	b = (x-x_[klo])/h; 
	ier_ = false;
	// Cubic spline polynomial is now evaluated.
	return a*y_[klo]+b*y_[khi]+((a*a*a-a)*y2_[klo]+(b*b*b-b)*y2_[khi])*(h*h)/6.0;
    };

    // Given arrays x[1..n] and y[1..n] containing a tabulated function, 
    // i.e., yi = f(xi), with x1 < x2 < .. . < xN, and given values yp1 
    // and ypn for the first derivative of the interpolating function at 
    // points 1 and n, respectively, this routine returns an array y2[1..n] 
    // that contains the second derivatives of the interpolating function at 
    // the tabulated points xi. If yp1 and/or ypn are equal to 1  1030 or 
    // larger, the routine is signaled to set the corresponding boundary
    // condition for a natural spline, with zero second derivative on that 
    // boundary.
    void spline(
        T yp1 = 1e30, 
        T ypn = 1e30
    ) {
        int i,k;
        T p,qn,sig,un,*u;
        if ( yp1 > 0.99e30 ) 
	    // The lower boundary condition is set either to be natural
            y2_[0] = u_[0] = 0.0;
        else { 
	    // or else to have a specified first derivative.
            y2_[0] = -0.5;
            u_[0] = (3.0/(x_[1]-x_[0]))*((y_[1]-y_[0])/(x_[1]-x_[0])-yp1);
        }
        // This is the decomposition loop of the tridiagonal algorithm.  
        // y2 and u are used for temporary storage of the decomposed factors.
        for ( i = 1; i < n_-1; i++ ) { 
	    sig = (x_[i]-x_[i-1])/(x_[i+1]-x_[i-1]);
	    p = sig*y2_[i-1]+2.0;
	    y2_[i] = (sig-1.0)/p;
	    u_[i] = (y_[i+1]-y_[i])/(x_[i+1]-x_[i]) - (y_[i]-y_[i-1])/(x_[i]-x_[i-1]);
	    u_[i] = (6.0*u_[i]/(x_[i+1]-x_[i-1])-sig*u_[i-1])/p;
        }
        if ( ypn > 0.99e30 ) 
	    // The upper boundary condition is set either to be natural 
	    qn = un = 0.0;
        else { 
	    // or else to have a specified first derivative.
	    qn = 0.5;
	    un = (3.0/(x_[n_]-x_[n_-1]))*(ypn-(y_[n_]-y_[n_-1])/(x_[n_]-x_[n_-1]));
        }
        y2_[n_] = (un-qn*u_[n_-1])/(qn*y2_[n_-1]+1.0);

        // This is the backsubstitution loop of the tridiagonal algorithm. 
        for ( k = n_-2; k >= 0; k-- ) 
        {
	    y2_[k] = y2_[k]*y2_[k+1]+u_[k];
        }
    };

    virtual void setValues( T x[], T y[], int n) 
    { 
	if ( n > n_ ) 
	{
	    delete[] y2_;
	    delete[] u_;
	    y2_ = new T[n+1];
	    u_ = new T[n+1];
	}
	Interpolator<T>::setValues(x,y,n);
	spline();
    };

  protected:

    // not allowed
    SplineInterpolator(SplineInterpolator<T>&);
    SplineInterpolator& operator=(SplineInterpolator<T>&);

  protected:

    T* y2_;
    T* u_;
};

}

#endif // _INTERPOLATOR_H
