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

#ifndef _OPD_XCOSM_H
#define _OPD_XCOSM_H

#include "opdBase.h"
#include <iostream>

namespace cosm {

template<typename T>
class opdXcosm : public opdBase<T> {

  public:

    opdXcosm(
        T ts_in,	// specimen thickness
        T tid_in,	// immersion thickness design
        T tia_in,	// immersion thickness actual  (not used)
        T tgd_in,	// coverglass thickness design
        T tga_in,	// coverglass thickness actual
        T ns_in,	// specimen refractive index 
        T nid_in,	// immersion refractive index design
        T nia_in,	// immersion refractive index actual
        T ngd_in,	// coverglass refractive index design
        T nga_in,	// coverglass refractive index actual
	T otd_in,	// optical tube length design
	T ota_in,	// optical tube length actual
	T mt_in,	// laterl magification
	T na_in		// numerical aperture
    ) : opdBase<T>(ts_in,((nia_in==nid_in)?0.0:tid_in),tia_in,tgd_in,tga_in,ns_in,nid_in,nia_in,ngd_in,nga_in), 
	t1_((otd_in-ota_in)*otd_in/(mt_in*mt_in-na_in*na_in)*ota_in), t2_(t1_/(2.0*nia_in)),
	niaons_(nia_in/ns_in), niaonid_(nia_in/nid_in), 
	niaongd_(nia_in/ngd_in), niaonga_(nia_in/nga_in),
	s4_((this->nssq_-this->niasq_)/this->ns_), 
	i4_((this->nidsq_-this->niasq_)/this->nid_), 
	g5_(tga_in*nga_in-tgd_in*ngd_in-this->niasq_*(tga_in/nga_in-tgd_in/ngd_in))
    {
	this->symmetric_ = this->symmetric_ && (otd_in == ota_in);
    };

    virtual ~opdXcosm() {};

    virtual T operator()( T rhoNAsq ) { 
	T tmp = sqrt(this->niasq_-rhoNAsq);
/*
	printf("ts: %e, term: %e, tmp: %e, s4: %e, nssq: %e, niaons: %e\n",
	        ts_, sqrt(nssq_-rhoNAsq), tmp, s4_, nssq_, niaons_);
	printf("phidz: %e, phis: %e, phitg: %e\n",
		(z_+t1_)*(tmp-nia_)+t2_*rhoNAsq,
	        ts_*(sqrt(nssq_-rhoNAsq)-niaons_*tmp-s4_),
	    -tid_*(sqrt(nidsq_-rhoNAsq)-niaonid_*tmp-i4_)
	    + tga_*(sqrt(ngasq_-rhoNAsq)-niaonga_*tmp)
	    -(tgd_*(sqrt(ngdsq_-rhoNAsq)-niaongd_*tmp)+g5_)
	   ); 
*/
	return (  
/*
	      (this->z_+t1_)*(tmp-this->nia_)+t2_*rhoNAsq
	    + this->ts_*(sqrt(this->nssq_-rhoNAsq)-niaons_*tmp-s4_)
	    - this->tid_*(sqrt(this->nidsq_-rhoNAsq)-niaonid_*tmp-i4_)
	    + this->tga_*(sqrt(this->ngasq_-rhoNAsq)-niaonga_*tmp)
	    -(this->tgd_*(sqrt(this->ngdsq_-rhoNAsq)-niaongd_*tmp)+g5_)

  	          09/24/2011 Daqi Dong
		  Removing terms of this->nia_, s4_, i4_, and g5_
		  because of test
         */
	      (this->z_+t1_)*tmp+t2_*rhoNAsq
	    + this->ts_*(sqrt(this->nssq_-rhoNAsq)-niaons_*tmp)
	    - this->tid_*(sqrt(this->nidsq_-rhoNAsq)-niaonid_*tmp)
	    + this->tga_*(sqrt(this->ngasq_-rhoNAsq)-niaonga_*tmp)
	    -(this->tgd_*(sqrt(this->ngdsq_-rhoNAsq)-niaongd_*tmp))

	);
    };
  
    T amplitude( T rhoNAsq ) {
	if ( this->niasq_ < rhoNAsq || this->nssq_ < rhoNAsq || 
	     this->nidsq_ < rhoNAsq || this->ngasq_ < rhoNAsq || this->ngdsq_ < rhoNAsq ) 
	{
	    return 0.0; 
	}
	return 1.0;
    };

    void dump() {
	std::cout.precision(12);
	std::cout.width(12);
	std::cout <<"S1: "<< this->ts_ <<", S2: "<< this->nssq_ <<", S3: "<< niaons_ <<", S4: "<< s4_ << std::endl;
	std::cout <<"G1: "<< niaongd_ <<", G2: "<< this->ngasq_ <<", G3: "<< this->ngdsq_ <<", G4: "<< niaonga_ <<", G5: "<< g5_ << std::endl;
	std::cout <<"I1: "<< this->tid_ <<", I2: "<< this->nidsq_ <<", I3: "<< niaonid_ <<", I4: "<< i4_ <<", I: "<< this->niasq_ << std::endl;
	std::cout <<"TUBE1: "<< t1_ <<", TUBE2: "<< t2_<< std::endl;
	std::cout <<"TT: "<< this->tga_ <<", TT0: "<< this->tgd_<<std::endl;
	std::cout <<"DZ: "<< this->z_<< ", DZ1: "<< this->nia_<< std::endl;
    };

  protected:

    // not allowed
    opdXcosm(opdXcosm<T>&);
    opdXcosm& operator=(opdXcosm<T>&);

  private:

    T t1_;
    T t2_;
    T niaons_;
    T niaonid_;
    T niaongd_;
    T niaonga_;
    T s4_;
    T i4_;
    T g5_;
};

}

#endif // _OPD_XCOSM_H
