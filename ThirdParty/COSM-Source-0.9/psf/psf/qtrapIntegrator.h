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
 * Returns the integral of the function func from a to b. The parameters EPS 
 * can be set to the desired fractional accuracy and JMAX so that 2 to the 
 * power JMAX-1 is the maximum allowed number of steps. Integration is 
 * performed by the trapezoidal rule.
 */
#ifndef _QTRAP_INTEGRATOR_H
#define _QTRAP_INTEGRATOR_H

#include "trapzdIntegrator.h"

namespace cosm {

template<typename T>
class QTrapIntegrator: public TrapzdIntegrator<T> {

  public:

    QTrapIntegrator( Functor<T>* func, T a, T b, T eps, int jmax ) 
	: TrapzdIntegrator<T>(func, a, b, eps, jmax) {};
    virtual ~QTrapIntegrator() {};

    // Returns integration value of a functor over interval [a,b]
    virtual T operator()() 
    {
        T  s;
        T olds = 0.0; // Initial value of olds is arbitrary. 
        for ( int j = 1; j <= this->jmax_; j++ ) 
        {
	    s = this->trapzd(j);
	    if (j > 5)  // Avoid spurious early convergence. 
	    {
	        if ( abs(s - olds) < this->eps_ * abs(olds) || 
		     (s == 0.0 && olds == 0.0) )
	        {
    		    this->ier_ = false;
		    return s;
	        }
	    }
	    olds = s;
        }
        this->ier_ = true;
        return 0.0; 
    };

  protected:

    // not allowed
    QTrapIntegrator(QTrapIntegrator<T>&);
    QTrapIntegrator& operator=(QTrapIntegrator<T>&);

};

}

#endif // _QTRAP_INTERGRATOR_H
