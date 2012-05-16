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

#include <iostream>

namespace cosm {

template<typename T>
void RadialPSF<T>::evaluate(
    PsfUser* user
)
{
    int Nz = psf_.length(0);
    int Nr = psf_.length(1);
    T zdist = 0;
    T rdist = 0;
    // the first half of z dimension has postive z, the second half negative
    int halfNz = Nz/2;
    // only need to calculate positive z if symmetric
    int zMax = symmetric_ ? halfNz : Nz-1;
    for ( int z = 0; z <= zMax; z++ ) 
    {
        zdist = (z > halfNz) ? (z-Nz) * deltaZ_ : z * deltaZ_;
        //std::cout <<"RadialPSF::evaluate; plane: "<<z<<", dz: "<<zdist<< std::endl;
        for ( int r = 0; r < Nr; r++ ) 
        {
	        rdist = r * deltaR_;
            complex<T> value = (*functor_)(zdist,rdist);
            psf_(z,r) = norm(value);
	        psfRe_(z,r) = value.real();
            psfIm_(z,r) = value.imag();
            if ( symmetric_ && z > 0 ) 
            {
                psf_(Nz-z, r) = psf_(z,r);
                psfRe_(Nz-z, r) = psfRe_(z,r);
                psfIm_(Nz-z, r) = 0 - psfIm_(z,r);
            }
        }

        //std::cout <<z<<" "<<psf_(z,Range::all())<< std::endl;
        //std::cout <<z+halfNz<<" "<<psf_(z,Range::all())<< std::endl;
        if ( user != NULL ) 
        {
            user->update(PsfUser::RADIAL, z, zMax);
	    }
	}
}

}
