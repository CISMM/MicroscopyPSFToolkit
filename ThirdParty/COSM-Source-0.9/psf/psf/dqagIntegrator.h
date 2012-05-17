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

#ifndef _DQAG_INTEGRATOR_H
#define _DQAG_INTEGRATOR_H

#include "integrator.h"
#include "quadpack/cquadpak.h"

namespace cosm {

template<typename T>
class DqagIntegrator : public Integrator<T> {

  public:

    DqagIntegrator( Functor<T>* func, T a, T b, T eps, int irule = 6)
	: Integrator<T>(func, a, b), eps_(eps), irule_(irule) {};
    virtual ~DqagIntegrator() {};


    // Returns integration value of a functor over interval [a,b]
    virtual T operator()()
    {
	int ier;
	double abserr;
	double value = 0;

	value = dqag(dqagFunction, this, double(this->a_),
                     double(this->b_), double(eps_), double(eps_), irule_,
                     &abserr, &neval_, &ier);
	error_ = T(abserr);
	ier_ = ( ier > 0 ) ? ier_ = true : ier_ = false;
	return T(value);
    };

    T errorValue() { return error_; };

  protected:

  static double dqagFunction( double val, void * cbData )
    {
      DqagIntegrator<T> * integrator = (DqagIntegrator<T> *)cbData;
      double value = (double)(*(integrator->func_))( (double) val );
      return value;
    };

  protected:
    // not allowed
    DqagIntegrator(DqagIntegrator<T>&);
    DqagIntegrator& operator=(DqagIntegrator<T>&);

  protected:

    T eps_;
    int irule_;
    int neval_;
    bool ier_;
    T error_;

};

};
#endif // _DQAG_INTEGRATOR_H
