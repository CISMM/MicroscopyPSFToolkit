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

#ifndef _BILINEAR_INTERPOLATOR_H
#define _BILINEAR_INTERPOLATOR_H

#include "interpolator.h"

namespace cosm {

template<typename T>
class BilinearInterpolator : public Interpolator<T> {

  public:

    BilinearInterpolator( T x[], T y[], int n ) 
	: Interpolator<T>(x, y, n) {};

    ~BilinearInterpolator() {};

    // Returns interpolation value
    virtual T operator()( T x ) { 
	this->ier_ = false;
	int ns = 0;
	int i = 0;
	for ( i = 0; i < this->n_; i++ )
	{	
	    if ( this->x_[i] == x )
	    {
	        return this->y_[i];
	    }
	    if ( this->x_[i] > x ) 
	    {
		ns = i - 1;
		break;
	    }
	}
	if ( ns < 0 || ns >= this->n_-1 )
	{
	    this->ier_ = true;
	    std::cout <<"x: "<<x<<", ns: "<<ns<<" x["<<i<<"]="<<this->x_[i]<<std::endl;
	    return 0;
	}
	return  (this->y_[ns+1] * (this->x_[ns+1]-x) + this->y_[ns] * (x-this->x_[ns]))/(this->x_[ns+1]-this->x_[ns]);
    };

    virtual void setValues( T x[], T y[], int n) { 
	Interpolator<T>::setValues(x, y, n); 
    };

  protected:

    // not allowed
    BilinearInterpolator(BilinearInterpolator<T>&);
    BilinearInterpolator& operator=(BilinearInterpolator<T>&);

  protected:

};

}

#endif // _BILINEAR_INTERPOLATOR_H
