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

#ifndef _ESTIMATE_LLS_H
#define _ESTIMATE_LLS_H

#include "est/estimateLinear.h"

namespace cosm {

template<typename T, int N>
class EstimateLLS : public EstimateLinear<T,N> {

  public:
    
    // Constructor
    EstimateLLS(
        Array<T,N>& psf,
        Array<T,N>& img,
	T pval
    ) : EstimateLinear<T,N>(psf, img, pval) { };

    // Destructor
    virtual ~EstimateLLS() {};

    // Calculate estimate
    virtual int run() 
    {
        this->estF_ = where(norm(this->psfF_) > this->val_, (this->imgF_*conj(this->psfF_))/norm(this->psfF_), 0);
	this->est_ = inverseFFT(this->estF_);	
        this->est_ = where( this->est_ < 0, 0, this->est_);
        return 0;
    };


    // Function to set pval
    void setPval( T pval ) { this->val_ = pval; };

    // Function to get pval
    T getPval() { return this->val_; };

  private:

    // not allowed
    EstimateLLS( EstimateLLS<T,N>& );
    void operator=( EstimateLLS<T,N>& );

};

}

#endif // _ESTIMATE_LLS_H
