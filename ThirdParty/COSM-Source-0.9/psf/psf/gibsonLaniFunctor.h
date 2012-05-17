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

#ifndef _GIBSON_LANI_FUNCTOR_H
#define _GIBSON_LANI_FUNCTOR_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "opdBase.h"
#include "bessel0nr.h"
#include "expFunctor.h"
#include <cmath>

namespace cosm {

template<typename T>
class GibsonLaniFunctor : public Functor<T> {

  public:

    GibsonLaniFunctor(T lambda, opdBase<T>& opd_in) 
    : exp_(&cos_), opd_(opd_in), k_(M_PI*2.0/lambda), r_(0), pre_(0) {};
    ~GibsonLaniFunctor() {};

    virtual T operator()( T x ) 
    { 
	T amplitude = opd_.amplitude(x);
	return amplitude == 0 ? 0 :
		j0_(pre_* sqrt(x))*(*exp_)(k_*opd_(x))*opd_.amplitude(x);
    };
    void setCos() { exp_ = &cos_; };
    void setSin() { exp_ = &sin_; };
    void setR( T r ) { r_ = r; pre_ = k_ * r_; };
    void setZ( T z ) { opd_.setZ(z); };
    T getK() { return k_; };
    T getR() { return r_; };
    T getPre() { return pre_; };
    T opd( T x ) { return opd_(x); };
  
  protected:

    // not allowed
    GibsonLaniFunctor( GibsonLaniFunctor<T>& );
    GibsonLaniFunctor& operator=( GibsonLaniFunctor<T>& );

  protected:
    
    J0nr<T> j0_;
    CosFunctor<T> cos_;
    SinFunctor<T> sin_;
    Functor<T>* exp_;
    opdBase<T>& opd_;
    T k_;
    T r_;
    T pre_; 

};

}

#endif // _GIBSON_LANI_FUNCTOR_H
