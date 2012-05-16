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

#ifndef _RADIAL_PSF_H
#define _RADIAL_PSF_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "psf/psfUser.h"
#include "psf/functorComplex.h"
#include "psf/bilinearInterpolator.h"
//#include "psf/splineInterpolator.h"
//#include "psf/ratInterpolator.h"
#include <blitz/array.h>

using namespace blitz;
namespace cosm {

enum EvalType {
    MAGNITUDE = 0x1,
    REAL      = 0x2,
    IMAGINARY = 0x4
};

   

template<typename T>
class RadialPSF {

  public:

    RadialPSF(
        int Nxy,
        int Nz,
        T deltaXY,
        T deltaZ,
        T deltaXYNyq,
        FunctorComplex<T>* functor,
        bool symmetric,
        int maxOversampling = 32
    ) : functor_(functor),
	psf_(Nz, Nxy), 
	psfRe_(Nz, Nxy),
	psfIm_(Nz, Nxy),
    nR_(0), nXY_(Nxy), nZ_(Nz),
	deltaR_(0), deltaXY_(deltaXY), deltaZ_(deltaZ), 
	deltaXYNyq_(deltaXYNyq),
	maxOversampling_(maxOversampling), 
	oversampling_(1), undersampled_(false), symmetric_(symmetric),
	interpolator_(0,0,0), z_(-1)
    {
        sampling();
    };

    ~RadialPSF() {};

    void evaluate(
        PsfUser* user = NULL
    );

    Array<T,2> psf( void ) { return psf_; };
	Array<T,2> psfRe( void ) { return psfRe_; };
	Array<T,2> psfIm( void ) { return psfIm_; };
    int nZ( void ) { return nZ_; };
    int nR( void ) { return nR_; };
    int nXY( void ) { return nXY_; };
    T deltaR( void ) { return deltaR_; };
    T deltaZ( void ) { return deltaZ_; };
    T deltaXY( void ) { return deltaXY_; };
    T deltaXYNyq( void ) { return deltaXYNyq_; };
    bool isSymmetric( void ) { return symmetric_; };
    bool isUndersampled( void ) { return undersampled_; };
    int oversampling( void ) { return oversampling_; };
    int maxOversampling( void ) { return maxOversampling_; };

    T exactValue( int z, int y, int x, EvalType type = MAGNITUDE )
    {
        T sq = deltaXY_ * deltaXY_;
        T r = sqrt(y*y*sq+x*x*sq);
        T zdist = (z > nZ_/2) ? (z-nZ_)*deltaZ_ : z*deltaZ_;
		complex<T> value = (*functor_)(zdist, r);
		switch ( type )
		{
		    case MAGNITUDE: return norm(value);
			case REAL: return real(value);
			case IMAGINARY: return imag(value);
		}
		return 0;
    };
                                                                                
    T interpolatedValue( int z, int y, int x, int oversampling = 1, EvalType type = MAGNITUDE )
    {
        T sq = (deltaXY_ * deltaXY_)/T(oversampling * oversampling);
        T r = sqrt(y*y*sq+x*x*sq);
        if ( z != z_ ) {
            z_ = z;
            slice_.resize(nR_);
			switch (type) 
			{
			    case MAGNITUDE: slice_ = psf_(z,Range::all()); break;
			    case REAL: slice_ = psfRe_(z,Range::all()); break;
			    case IMAGINARY: slice_ = psfIm_(z,Range::all()); break;
            }
            T* ra = new T[nR_];
            for ( int i = 0; i < nR_; i++ )
            {
                ra[i] = i*deltaR_;
            }
            interpolator_.setValues((T*)(ra), (T*)(slice_.data()), nR_);
        }
	T val = interpolator_(r);
        return val;
    };

  private:

    void sampling( void ) 
    { 
        // Oversampling is needed because of radial sweep. 1/4 with respect 
        // to the Nyquist distance seems to work well with the cubic spline 
        // interpolation that we use in the PSF computation. */
        deltaR_ =  deltaXYNyq_ / 4;

		/* Find how much we have to oversample by for the radial values */
        if (deltaXY_/deltaXYNyq_ > 1.0)
        {
            undersampled_ = true;
        }
        oversampling_ = (int)(deltaXY_/deltaR_ + 0.5);
 
		if ( oversampling_ <= 0 ) 
        {
	        oversampling_ = 1;
        }
        if ( oversampling_ > maxOversampling_)
        {
            oversampling_ = maxOversampling_;
		}

        // make oversampling power of 2 to match xcosm
        int i = 0;
        for ( i = 1; i < oversampling_; i = i << 1 );
        oversampling_ = i;

        deltaR_ = deltaXY_/oversampling_;

        // adjust the radial samples
        nR_ = (int)(((psf_.length(1)+1)/2 + 1.0)*M_SQRT2 + 1.0)*oversampling_ + 5;
        psf_.resize(psf_.length(0), nR_);
		psfRe_.resize(psf_.length(0), nR_);
		psfIm_.resize(psf_.length(0), nR_);
        cout <<"undersampled: "<< (undersampled_ ? "true" : "false") << endl;
        cout <<"oversampling: "<< oversampling_ << endl;
        cout <<"deltar: "<< deltaR_ <<", Nr: "<< nR_ << endl;
    };

  private:

    // not allowed
    RadialPSF( RadialPSF<T>& );
    RadialPSF& operator=( RadialPSF<T>& );

  private:

    FunctorComplex<T>* functor_;
    Array<T,2> psf_;
	Array<T,2> psfRe_;
	Array<T,2> psfIm_;
    int nR_;
    int nXY_;
    int nZ_;
    T deltaR_;
    T deltaXY_;
    T deltaZ_;
    T deltaXYNyq_;
    int maxOversampling_;
    int oversampling_;
    bool undersampled_;
    bool symmetric_;
    BilinearInterpolator<T> interpolator_;
    //SplineInterpolator<T> interpolator_;
    //RatInterpolator<T> interpolator_;
    int z_;
    Array<T,1> slice_;

};

}

#include "psf/radialPSF.cxx"

#endif  // _RADIAL_PSF_H
