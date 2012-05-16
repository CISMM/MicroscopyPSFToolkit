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

#ifndef _BESSEL_N_M_H
#define _BESSEL_N_M_H

#include "functor.h"
#include <cmath>

namespace cosm {

template<typename T>
class Jnm : public Functor<T> {

  public: 

    Jnm( int n_in = 0 ) : n_(n_in) {};
    ~Jnm() {};

    // Returns the Bessel function J0(x) for any x.
    T operator()( T x ) 
    {
	switch(n_) 
	{
	   case 0: return T(j0(x));
	   case 1: return T(j1(x));
	   default : return T(jn(n_,x));
	}
    }
    void n( int n_in ) { n_ = n_in; };

  protected:

    // not allowed
    Jnm(Jnm<T>&);
    Jnm& operator=(Jnm<T>&);

  protected:

    int n_;
};

}
#endif // _BESSEL_N_M_H
