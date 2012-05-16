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

#include "psf/qtrapIntegrator.h"
#include "psf/qsimpIntegrator.h"
#include "psf/qrombIntegrator.h"
#include "psf/dqagIntegrator.h"
#include <blitz/timer.h>
#include <iostream> 

using namespace cosm;

//#define REAL double
#define REAL float

template<typename T>
class TestFunctor : public Functor<T> {
  public:
    TestFunctor() : Functor<T>() {};
    ~TestFunctor() {};
    T operator()( T x ) { return cos(x)*j0(sqrt(x)); };
};

int main ( int argc, char* argv[] ) {

    REAL eps = 1e-5;
    int jmax = 20;

    const int max_int = 4;
    TestFunctor<REAL> functor;
    QTrapIntegrator<REAL> qtrap(&functor, 0.0, 1.0, eps, jmax);
    QSimpIntegrator<REAL> qsimp(&functor, 0.0, 1.0, eps, jmax);
    QRombIntegrator<REAL> qromb(&functor, 0.0, 1.0, eps, jmax);
    DqagIntegrator<REAL>  dqag(&functor, 0.0, 1.0, eps);
    Integrator<REAL>* integrator[max_int] = {&qtrap, &qsimp, &qromb, &dqag};

    blitz::Timer timer;
    REAL value = 0.0;
    for ( int i = 0; i < max_int; i++ ) 
    {
	timer.start();
	for ( int j = 0; j < 10000; j++ ) 
	{
            value = (*integrator[i])();
	}
	timer.stop();
        std::cout << "value: "<<value<<", time for integration: "<<timer.elapsedSeconds()<< std::endl;
    }
    return 0;
}
