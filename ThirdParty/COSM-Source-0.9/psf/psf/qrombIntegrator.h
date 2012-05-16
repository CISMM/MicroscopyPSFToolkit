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
 * Here EPS is the fractional accuracy desired, as determined by the 
 * extrapolation error estimate; JMAX limits the total number of steps; 
 * K is the number of points used in the extrapolation.
 * Returns the integral of the function func from a to b. Integration is 
 * performed by Rombergs method of order 2K, where, e.g., K=2 is Simpsons rule.
 */

#ifndef _QROMB_INTEGRATOR_H
#define _QROMB_INTEGRATOR_H

#include "trapzdIntegrator.h"
#include "polyInterpolator.h"
#include <iostream>

namespace cosm {

template<typename T>
class QRombIntegrator: public TrapzdIntegrator<T> {

  public:

    QRombIntegrator( Functor<T>* func, T a, T b, T eps, int jmax, int k = 5 ) 
	: TrapzdIntegrator<T>(func, a, b, eps, jmax), polint_(0,0,k), k_(k)
    {
        // These store the successive trapezoidal approximations and their 
        // relative stepsizes. 
        s_ = new T[jmax];
	h_ = new T[jmax+1]; 
    }
    virtual ~QRombIntegrator() 
    {
	delete[] s_;
	delete[] h_;
    };

    // Returns integration value of a functor over interval [a,b]
    virtual T operator()() 
    {
        T ss, dss;
        h_[0] = 1.0;
        for ( int j = 0; j < this->jmax_; j++ ) {
            s_[j] = this->trapzd(j);
            if ( j >= k_-1 ) {
                polint_.setValues(&h_[j-k_],&s_[j-k_],k_);
                ss = polint_(0.0);
		if ( polint_.error() ) {
		    std::cout <<"qrombIntegrator; error in polint"<<std::endl;
		    this->ier_ = true; 
		    return -1.0;
		}
		dss = polint_.errorValue();
                if ( abs(dss) <= this->eps_*abs(ss) ) 
	        {
		    this->ier_ = false;
		    return ss;
	        }
            }
            /* 
	     * This is a key step: The factor is 0.25 even though the stepsize 
	     * is decreased by only 0.5. This makes the extrapolation a 
	     * polynomial in h2 as allowed by equation (4.2.1), not just 
	     * a polynomial in h.
	     */
            h_[j+1] = 0.25 * h_[j];
        }
	std::cout <<"qrombIntegrator; error in convergence"<<std::endl;
        this->ier_ = true;
        return -2.0; 
    };

  protected:

    // not allowed
    QRombIntegrator(QRombIntegrator<T>&);
    QRombIntegrator& operator=(QRombIntegrator<T>&);

  protected:

    PolyInterpolator<T> polint_;
    int k_;
    T* s_;
    T* h_;

};

}

#endif // _QROMB_INTERGRATOR_H
