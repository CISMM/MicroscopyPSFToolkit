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
 * performed by Simpsons rule.
 */
#ifndef _QSIMP_INTEGRATOR_H
#define _QSIMP_INTEGRATOR_H

#include "trapzdIntegrator.h"

namespace cosm {

template<typename T>
class QSimpIntegrator: public TrapzdIntegrator<T> {

  public:

    QSimpIntegrator( Functor<T>* func, T a, T b, T eps, int jmax = 10 ) 
	: TrapzdIntegrator<T>(func, a, b, eps, jmax) {};
    virtual ~QSimpIntegrator() {};

    // Returns integration value of a functor over interval [a,b]
    virtual T operator()() 
    {
        T s, st, ost = 0.0, os = 0.0;

        for ( int j = 1; j <= this->jmax_; j++ ) 
	{
            st = this->trapzd(j);
            s = (4.0*st-ost)/3.0;    
            if ( j > 5 )    // Avoid spurious early convergence. */
	    {
                if ( abs(s-os) < this->eps_ * abs(os) || (s == 0.0 && os == 0.0) ) 
	        {
    		    this->ier_ = false;
		    return s;
	        }
	    }
            os = s;
            ost = st;
        }
        this->ier_ = true;
        return 0.0; 
    };

  protected:

    // not allowed
    QSimpIntegrator(QSimpIntegrator<T>&);
    QSimpIntegrator& operator=(QSimpIntegrator<T>&);

};

}

#endif // _QSIMP_INTERGRATOR_H
