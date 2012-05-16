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

#ifndef _GIBSON_LANI_PSF_FUNCTOR_H
#define _GIBSON_LANI_PSF_FUNCTOR_H

#include "psf/psfFunctor.h"
#include "psf/dqagIntegrator.h"
//#include "psf/qsimpIntegrator.h"
#include "psf/gibsonLaniFunctor.h"
#include "psf/opdXcosm.h"
#include <complex>

namespace cosm {

template<typename T>
class GibsonLaniPsfFunctor : public PsfFunctor<T> {

  public:

    GibsonLaniPsfFunctor(
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
        gibsonLaniFunctor_(lambda, opd_),
        integrator_(&gibsonLaniFunctor_, 0, na*na, absError)
    {};

    ~GibsonLaniPsfFunctor() {};

    virtual complex<T> operator()( T z, T r ) 
    {
        gibsonLaniFunctor_.setZ(z);
        gibsonLaniFunctor_.setR(r);
		gibsonLaniFunctor_.setCos();
        T realPart = integrator_();
        gibsonLaniFunctor_.setSin();
        T imagPart = integrator_();
        // cout <<"("<<z<<","<<r<<"), real: "<<realPart<<", imag: "<<imagPart<<endl;
        complex<T> retVal = std::complex<T>(realPart, imagPart);
        // std::cout <<"psf("<<z<<","<<r<<")="<<retVal<<std::endl;

	return retVal;
    };
   
    virtual bool isSymmetric() { return opd_.isSymmetric(); };
  
  protected:

    // not allowed
    GibsonLaniPsfFunctor( GibsonLaniPsfFunctor<T>& );
    GibsonLaniPsfFunctor& operator=( GibsonLaniPsfFunctor<T>& );

  protected:
    
    opdXcosm<T> opd_;
    GibsonLaniFunctor<T> gibsonLaniFunctor_;
    DqagIntegrator<T> integrator_;
    //QSimpIntegrator<T> integrator_;

};

}

#endif // _GIBSON_LANI_PSF_FUNCTOR_H
