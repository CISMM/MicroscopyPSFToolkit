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

#ifndef _HAEBERLE_PSF_FUNCTOR_H
#define _HAEBERLE_PSF_FUNCTOR_H

#include "psf/psfFunctor.h"
#include "psf/dqagIntegrator.h"
#include "psf/haeberleFunctor.h"
#include "psf/opdXcosm.h"
#include <complex>

namespace cosm {

template<typename T>
class HaeberlePsfFunctor : public PsfFunctor<T> {

  public:

    HaeberlePsfFunctor(
        T ts,           // specimen thickness
        T tid,          // immersion thickness design
        T tia,          // immersion thickness actual
        T tgd,          // coverglass thickness design
        T tga,          // coverglass thickness actual
        T ns,           // specimen refractive index
        T nid,          // immersion refractive index design
        T nia,          // immersion refractive index actual
        T ngd,          // coverglass refractive index design
        T nga,          // coverglass refractive index actual
        T tld,          // tube length design
        T tla,          // tube length actual
        T lm,           // lateral magnification
        T na,           // numerical aperture (NA)
        T lambda,
        T absError
    ) : PsfFunctor<T>(), 
	opd_(ts, tid, tia, tgd, tga, ns, nid, nia, ngd, nga, tld, tla, lm, na),
        haeberleFunctor_(lambda, opd_),
        integrator_(&haeberleFunctor_, 0, na*na, absError)
    {};

    ~HaeberlePsfFunctor() {};

    virtual complex<T> operator()( T z, T r ) 
    {
        complex<T> retVal = 0;
		T realPart = 0;
        T imagPart = 0;
        haeberleFunctor_.setZ(z);
        haeberleFunctor_.setR(r);

        haeberleFunctor_.setType(HaeberleFunctor<T>::HAEBERLE_I0);
        haeberleFunctor_.setCos();
        realPart = integrator_();
        haeberleFunctor_.setSin();
        imagPart = integrator_();
        retVal += std::complex<T>(realPart, imagPart);

        haeberleFunctor_.setType(HaeberleFunctor<T>::HAEBERLE_I1);
        haeberleFunctor_.setCos();
        realPart = integrator_();
        haeberleFunctor_.setSin();
        imagPart = integrator_();
        retVal += std::complex<T>(2 *realPart, 2 *imagPart);

        haeberleFunctor_.setType(HaeberleFunctor<T>::HAEBERLE_I2);
        haeberleFunctor_.setCos();
        realPart = integrator_();
        haeberleFunctor_.setSin();
        imagPart = integrator_();
        retVal += std::complex<T>(realPart, imagPart);

        // std::cout <<"psf("<<z<<","<<r<<")="<<retVal<<std::endl;
        return retVal;
    };
   
    virtual bool isSymmetric() { return opd_.isSymmetric(); };
  
  protected:

    // not allowed
    HaeberlePsfFunctor( HaeberlePsfFunctor<T>& );
    HaeberlePsfFunctor& operator=( HaeberlePsfFunctor<T>& );

  protected:
    
    opdXcosm<T> opd_;
    HaeberleFunctor<T> haeberleFunctor_;
    DqagIntegrator<T> integrator_;

};

}

#endif // _HAEBERLE_PSF_FUNCTOR_H
