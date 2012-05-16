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

#ifndef _INTERPOLATOR_H
#define _INTERPOLATOR_H

#include "functor.h"

namespace cosm {

template<typename T>
class Interpolator : public Functor<T> {

  public:

    Interpolator( T x[], T y[], int n ) 
	: Functor<T>(), x_(x), y_(y), n_(n), error_(0), ier_(false)  {};

    ~Interpolator() {};

    // Returns interpolation value
    virtual T operator()( T x ) = 0;

    virtual void setValues( T x[], T y[], int n) { x_ = x; y_ = y; n_ = n; };

    T errorValue() { return error_; };
    bool error() { return ier_; }; 

  protected:

    // not allowed
    Interpolator(Interpolator<T>&);
    Interpolator& operator=(Interpolator<T>&);

  protected:

    T* x_;
    T* y_;
    int n_;
    T error_;
    bool ier_;

};

}

#endif // _INTERPOLATOR_H
