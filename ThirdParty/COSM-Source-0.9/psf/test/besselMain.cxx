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

#include "psf/bessel0m.h"
#include "psf/bessel0nr.h"
#include "psf/bessel0nl.h"

#include <blitz/timer.h>
#include <iostream> 

#define _USE_MATH_DEFINES

using namespace cosm;

int main ( int argc, char* argv[] ) {

	int N = 100; 
	double x = 0.0;
	double xMax = 100.0;
	double delta = xMax/double(N);

	J0m<double> j0m;
	J0nr<double> j0nr;
	J0nl<double> j0nl;

	std::cout.precision(12);
	std::cout.width(12);

	while ( x <= xMax ) {
	    x += delta;
	    std::cout <<"x: "<< x <<", m: " << j0(x) <<", nr: "<< j0nr(x) <<", nl: "<< j0nl(x) << std::endl;
	}

	x = 0.0;
	xMax = 10;
	delta = xMax/1000000;
	blitz::Timer timer;
	long double time;

	timer.start();
	while ( x <= xMax ) {
		x += delta;
        double y = 0;
		y = j0(x);
	}
	timer.stop();
	time = timer.elapsedSeconds();
	std::cout << "time for j0: "<< time << std::endl;

	x = 0.0;
    timer.start();
	while ( x <= xMax )
	{
		x += delta;
	    double y = 0;
		y = j0nr(x);
	}
	timer.stop();
    time = timer.elapsedSeconds();
	std::cout << "time for j0nr: "<< time << std::endl;

	x = 0.0;
	timer.start();
	while ( x <= xMax )
	{
	    x += delta;
	    double y = 0;
		y = j0nl(x);
	}
	timer.stop();
	time = timer.elapsedSeconds();
	std::cout << "time for j0nl: "<< time << std::endl;

	return 0;
}
