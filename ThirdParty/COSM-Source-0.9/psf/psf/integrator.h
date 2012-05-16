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

#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#include "functor.h"

namespace cosm {

template<typename T>
class Integrator {

  public: 

    Integrator( Functor<T>* func, T a, T b ) 
	: func_(func), a_(a), b_(b), ier_(false) {};
    virtual ~Integrator() {};

    // Returns integration value of a functor over interval [a,b]
    virtual T operator()() = 0;

    void integrand( Functor<T>* func ) { func_ = func; };
    void lower( T a ) { a_ = a; };
    void upper( T b ) { b_ = b; };
    bool error() { return ier_; };

  protected:

    // not allowed
    Integrator(Integrator<T>&);
    Integrator& operator=(Integrator<T>&);

  protected:

    Functor<T>* func_;
    T a_;
    T b_;
    bool ier_;

};

}

#endif // _INTEGRATOR_H
