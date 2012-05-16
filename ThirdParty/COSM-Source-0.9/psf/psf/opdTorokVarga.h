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

#ifndef _OPD_TOROK_VARGA_H
#define _OPD_TOROK_VARGA_H

#include "opdBase.h"

namespace cosm {

// Base class to compute the OPD for given optical parameters and 
// model by Torok and Varga, as presented in paper by O. Haeberle.
template<typename T>
class opdTorokVarga : public opdBase<T> {

  public:

    // constructor
    opdTorokVarga(
        T ts, 		// specimen thickness
        T tid,		// immersion thickness design
        T tia,		// immersion thickness actual
        T tgd,		// coverglass thickness design
        T tga,		// coverglass thickness actual
        T ns,		// specimen refractive index 
        T nid,		// immersion refractive index design
        T nia,		// immersion refractive index actual
        T ngd,		// coverglass refractive index design
        T nga		// coverglass refractive index actual
    ) : opdBase<T>(ts,tid,tia,tgd,tga,ns,nid,nia,ngd,nga), 
        nsts_(ns*ts), niatia_(nia*tia), 
	ngatga_(nga*tga), ngdtgd_(ngd*tgd),
	niaonssq_(nia*nia/ns/ns), niaonidsq_(nia*nia/nid/nid), 
	niaongdsq_(nia*nia/ngd/ngd), niaongasq_(nia*nia/nga/nga)
    {};

    // destructor
    ~opdTorokVarga() {};

    // function to evaluate the OPD
    virtual T operator()( T rhoNAsq ) { 
	return (  
		this->niasq_*term1(this->niasq_,rhoNAsq)  
	      +	ngatga_*(term1(this->ngasq_,rhoNAsq)-term(niaongasq_,this->niasq_,rhoNAsq)) 
	      -	ngdtgd_*(term1(this->ngdsq_,rhoNAsq)-term(niaongdsq_,this->niasq_,rhoNAsq)) 
	      -	nidtid_*(term1(this->nidsq_,rhoNAsq)-term(niaonidsq_,this->niasq_,rhoNAsq)) 
	      -	nsts_  *(term1(this->nssq_,rhoNAsq)-term(niaonssq_, this->niasq_,rhoNAsq)) 
	);
    };

  protected:

    // not allowed
    opdTorokVarga(opdTorokVarga<T>&);
    opdTorokVarga& operator=(opdTorokVarga<T>&);

  private:

    // pre-computed constants

    T nsts_;	// ns * ts
    T niatia_;	// nia * tia
    T nidtid_;	// nid * tid
    T ngatga_;  // nga * tga
    T ngdtgd_;  // ngd * tgd

    T niaonssq_;  // (nia/ns)**2
    T niaonidsq_; // (nia/nid)**2
    T niaongdsq_; // (nia/ngd)**2
    T niaongasq_; // (nia/nga)**2

};

}

#endif // _OPD_TOROK_VARGA_H
