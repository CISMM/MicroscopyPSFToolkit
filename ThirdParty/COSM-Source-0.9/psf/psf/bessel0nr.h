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

// Functor that returns the Bessel function J0(x) for any x.
// Based on code in Numerical Recipies for C.

#ifndef _BESSEL_0_NR_H
#define _BESSEL_0_NR_H

#include "functor.h"
#include <math.h>

namespace cosm {

template<typename T>
class J0nr : public Functor<T> {

  public:

    J0nr() {};
    virtual ~J0nr() {};

    // Returns the Bessel function J0(x) for any x.
    T operator()( T x ) 
    {
        // Accumulate polynomials in double precision.
        T xx;
        T ans;
	T xabs = x < 0 ? -x : x;
        if ( xabs < 8.0 ) {  
	    /* Direct rational function fit. */
            T y = x*x;
            T ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7
    	           + y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
            T ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718
	           + y*(59272.64853+y*(267.8532712+y*1.0))));
            ans = ans1/ans2;
        } else {  
	    /* Fitting function (6.5.9). */
	    T z = 8.0/xabs;
	    T y = z*z;
	    xx = xabs-0.785398164;
	    T ans1 = 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
	           + y*(-0.2073370639e-5+y*0.2093887211e-6)));
	    T ans2 = -0.1562499995e-1+y*(0.1430488765e-3
	           + y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934945152e-7)));
	    ans = sqrt(0.636619772/xabs)*(cos(xx)*ans1-z*sin(xx)*ans2);
        }
        return ans;
    }

  protected:

    // not allowed
    J0nr(J0nr<T>&);
    J0nr& operator=(J0nr<T>&);

};

}

#endif // _BESSEL_0_NR_H
