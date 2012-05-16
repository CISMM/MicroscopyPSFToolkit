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

#ifndef _OPD_GIBSON_LANI_H
#define _OPD_GIBSON_LANI_H

#include "opdBase.h"

namespace cosm {

template<typename T>
class opdGibsonLani : public opdBase<T> {

  public:

    opdGibsonLani(
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
    ) : opdBase<T>(ts,tid,tia,tgd,tga,ns,nid,nia,ngd,nga) 
    {};

    ~opdGibsonLani() {};

    virtual T operator()( T rhoNAsq ) { 
	return (  term(this->tgd_,this->ngdsq_,rhoNAsq)  
		- term(this->tga_,this->ngasq_,rhoNAsq)  
		+ term(this->tid_,this->nidsq_,rhoNAsq)  
		- term(this->tia_,this->niasq_,rhoNAsq)  
		+ term(this->ts_, this->nssq_, rhoNAsq) );
    };

  protected: 

    // not allowed
    opdGibsonLani(opdGibsonLani<T>&);
    opdGibsonLani& operator=(opdGibsonLani<T>&);

};

}

#endif // _OPD_GIBSON_LANI_H
