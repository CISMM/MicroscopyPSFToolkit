/****************************************************************************
 * Copyright (c) 2004 Einir Valdimarsson and Chrysanthe Preza
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

#include "psf/gibsonLaniPsfFunctor.h"
#include "psf/bilinearInterpolator.h"
#include "psf/polyInterpolator.h"
#include "psf/ratInterpolator.h"
#include "psf/splineInterpolator.h"
#include <blitz/timer.h>
#include <iostream> 

using namespace cosm;

#define REAL float

template<typename T>
class TestFunctor : public Functor<T> {
  public:
    TestFunctor() : Functor<T>() {};
    ~TestFunctor() {};
    //T operator()( T x ) { return cos(x)*j0(sqrt(x)); };
    T operator()( T x ) { return cos(x); };
};

int main ( int argc, char* argv[] ) {

    const int n = 100;
    const int k = 2;
    REAL x[n+1];
    REAL y[n+1];
    //const int max_int = 4;
    const int max_int = 3;
    REAL z[max_int][k*(n+1)];
    //TestFunctor<REAL> functor;
    REAL lambda = 0.000580;
    GibsonLaniPsfFunctor<REAL> functor( 
        0.050,          // specimen thickness
        0.16,           // immersion thickness design
        0.16,           // immersion thickness actual
        0.170,          // coverglass thickness design
        0.120,          // coverglass thickness actual
        1.33,           // specimen refractive index
        1.0,            // immersion refractive index design
        1.0,            // immersion refractive index actual
        1.525,          // coverglass refractive index design
        1.525,          // coverglass refractive index actual
	0.160,		// otld	
	0.160,		// otla	
	1.0,		// mt
        0.9,            // numerical aperture
	lambda,
	1e-6
    );
    for ( int i = 0; i <= n; i++ ) 
    {
	x[i] = i/REAL(n);
	//y[i] = functor(x[i]);
	y[i] = functor(0,x[i]);
    }
    BilinearInterpolator<REAL> bilinear(x,y,n+1);
    PolyInterpolator<REAL> poly(x,y,n+1);
    RatInterpolator<REAL> rat(x,y,n+1);
    SplineInterpolator<REAL> spline(x,y,n+1);
    //Interpolator<REAL>* interpolator[max_int] = {&bilinear, &poly, &rat, &spline};
    Interpolator<REAL>* interpolator[max_int] = {&bilinear, &rat, &spline};

    blitz::Timer timer;
    REAL value = 0.0;
    for ( int i = 0; i < max_int; i++ ) 
    {
	timer.start();
	for ( int j = 0; j < k*n; j++ ) 
	{
            z[i][j] = (*interpolator[i])(j/REAL(k*n));
	    if ( interpolator[i]->error() )
	    {
		std::cout <<"Error "<<i<<std::endl;
	    }
	}
	timer.stop();
        std::cout <<i<< ", value: "<<value<<", time for interpolation: "<<timer.elapsedSeconds()<< std::endl;
    }
    for ( int j = 0; j < k*n; j++ ) 
    {
	//std::cout << functor(j/REAL(k*n)) <<" "; 
	std::cout << functor(0,j/REAL(k*n)) <<" "; 
        for ( int i = 0; i < max_int; i++ ) 
	{
	    std::cout << z[i][j] <<" "; 
	}
	
	std::cout << std::endl;
    }

    return 0;
}
