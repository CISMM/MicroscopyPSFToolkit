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
 * This routine computes the nth stage of refinement of an extended trapezoidal
 * rule. func is input as a pointer to the function to be integrated between 
 * limits a and b, also input. When called with n=1, the routine returns the 
 * crudest estimate of int_b_a f(x)dx. Subsequent calls with n=2,3,...
 * (in that sequential order) will improve the accuracy by adding 2n-2 
 * additional interior points.
 */
#ifndef _TRAPZD_INTEGRATOR_H
#define _TRAPZD_INTEGRATOR_H

#include "integrator.h"

namespace cosm {

template<typename T>
class TrapzdIntegrator: public Integrator<T> {

  public:

    TrapzdIntegrator( Functor<T>* func, T a, T b, T eps, int jmax ) 
    : Integrator<T>(func,a,b), s_(0), eps_(eps), jmax_(jmax) {};
    virtual ~TrapzdIntegrator() {};

    // Returns integration value of a functor over interval [a,b]
    virtual T operator()() = 0;

  protected: 

    T trapzd( int n )
    {
        int it, j;
        if ( n == 1 ) {
	    s_ = 0.5 * (this->b_ - this->a_) * ((*this->func_)(this->a_) + (*this->func_)(this->b_));
        } else {
	    for ( it = 1, j = 1; j < n-1; j++ ) {
	        it <<= 1;
	    }
	    T tnm = T(it);
	    /* This is the spacing of the points to be added. */
	    T del = (this->b_ - this->a_)/tnm; 
	    T x = this->a_ + 0.5 * del;
	    T sum = 0.0;
	    for ( j = 0; j < it; j++, x += del ) {
	        sum += (*this->func_)(x);
	    }
	    /* This replaces s by its refined value. */
	    s_ = 0.5 * (s_ + (this->b_ - this->a_) * sum/tnm); 
        }
	return s_;
    }

    T abs(T x) { return x < 0 ? -x : x; };

  protected:

    // not allowed
    TrapzdIntegrator(TrapzdIntegrator<T>&);
    TrapzdIntegrator& operator=(TrapzdIntegrator<T>&);

  protected:

    T s_;
    T eps_;
    int jmax_;

};

}

#endif // _INTERGRATOR_H
