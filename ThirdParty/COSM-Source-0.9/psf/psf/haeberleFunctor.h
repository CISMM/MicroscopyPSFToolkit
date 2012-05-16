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

#ifndef _HAEBERLE_FUNCTOR_H
#define _HEABERLE_FUNCTOR_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "opdBase.h"
#include "besselnm.h"
#include "expFunctor.h"
#include <cmath>

namespace cosm {

template<typename T>
class HaeberleFunctor : public Functor<T> {

  public:

    enum Type {
	HAEBERLE_I0 = 0x0,
	HAEBERLE_I1 = 0x1,
	HAEBERLE_I2 = 0x2
    };

  public:

    HaeberleFunctor(T lambda, opdBase<T>& opd_in) 
    : exp_(&cos_), opd_(opd_in), k_(M_PI*2.0/lambda), r_(0), pre_(0),
      ns_(opd_.ns()), nia_(opd_.nia()), nga_(opd_.nga()), 
      ns2_(ns_*2), nia2_(nia_*2), nga2_(nga_*2),
      niasq_(nia_*nia_), niaons_(nia_/ns_), niaonga_(nia_/nga_), 
      type_(HAEBERLE_I0)
    { };
    ~HaeberleFunctor() {};

    virtual T operator()( T x ) 
    { 
	if ( opd_.amplitude(x) == 0 )   
 	{
	    return 0;
	}
	switch( type_ ) 
	{
	    case HAEBERLE_I0: 
		return integrand0(x);
	    case HAEBERLE_I1: 
		return integrand1(x);
	    case HAEBERLE_I2: 
		return integrand2(x);
	}
	return 0;
    };

    T integrand0( T x ) 
    {	
	jn_.n(0);
        T ctheta1 = sqrt(1 - x/niasq_);
        T stheta1 = sqrt(1 - ctheta1*ctheta1);
        T stheta2 = stheta1 * niaonga_;
        T ctheta2 = sqrt(1 - stheta2*stheta2);
        T stheta3 = stheta1 * niaons_;
        T ctheta3 = sqrt(1 - stheta3*stheta3);
                                                                                
        T retVal = nia2_*ctheta1/(nia_*ctheta1+nga_*ctheta2) *  
		   nga2_*ctheta2/(nga_*ctheta2+ns_*ctheta3);
        retVal +=  nia2_ *ctheta1/(nga_*ctheta1+nia_*ctheta2) *
	           nga2_*ctheta2/(ns_*ctheta2+nga_*ctheta3) * ctheta3;
        retVal *= jn_(pre_* sqrt(x))*(*exp_)(k_*opd_(x))*opd_.amplitude(x)*sqrt(ctheta1);
	return retVal;
    };

    T integrand1( T x ) 
    {
	jn_.n(1);
        T ctheta1 = sqrt(1 - x/niasq_);
        T stheta1 = sqrt(1 - ctheta1*ctheta1);
        T stheta2 = stheta1 * niaonga_;
        T ctheta2 = sqrt(1 - stheta2*stheta2);
        T stheta3 = stheta1 * niaons_;
        T ctheta3 = sqrt(1 - stheta3*stheta3);
                                                                                
        T retVal = nia2_ *ctheta1/(nga_*ctheta1+nia_*ctheta2) *
	           nga2_*ctheta2/(ns_*ctheta2+nga_*ctheta3) * stheta3;
        retVal *= jn_(pre_* sqrt(x))*(*exp_)(k_*opd_(x))*opd_.amplitude(x)*sqrt(ctheta1);
	return retVal;
    };

    T integrand2( T x ) 
    {
	jn_.n(2);
        T ctheta1 = sqrt(1 - x/niasq_);
        T stheta1 = sqrt(1 - ctheta1*ctheta1);
        T stheta2 = stheta1 * niaonga_;
        T ctheta2 = sqrt(1 - stheta2*stheta2);
        T stheta3 = stheta1 * niaons_;
        T ctheta3 = sqrt(1 - stheta3*stheta3);
                                                                                
        T retVal = nia2_*ctheta1/(nia_*ctheta1+nga_*ctheta2) *  
		   nga2_*ctheta2/(nga_*ctheta2+ns_*ctheta3);
        retVal -=  nia2_ *ctheta1/(nga_*ctheta1+nia_*ctheta2) *
	           nga2_*ctheta2/(ns_*ctheta2+nga_*ctheta3) * ctheta3;
        retVal *= jn_(pre_* sqrt(x))*(*exp_)(k_*opd_(x))*opd_.amplitude(x)*sqrt(ctheta1);
	return retVal;
    };

    void setType( Type type ) { type_ = type; };
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
    HaeberleFunctor( HaeberleFunctor<T>& );
    HaeberleFunctor& operator=( HaeberleFunctor<T>& );

  protected:
    
    Jnm<T> jn_;
    CosFunctor<T> cos_;
    SinFunctor<T> sin_;
    Functor<T>* exp_;
    opdBase<T>& opd_;
    T k_;
    T r_;
    T pre_; 
    T ns_;
    T nia_;
    T nga_;
    T ns2_;
    T nia2_;
    T nga2_;
    T niasq_;
    T niaons_;
    T niaonga_;
    Type type_;

};

}

#endif // _HAEBERLE_FUNCTOR_H
