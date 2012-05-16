
#include "psf/opdXcosm.h"
#include "psf/opdGibsonLani.h"
#include "psf/opdTorokVarga.h"
#include <iostream>

using namespace cosm;

int main(
    int argc,
    char* argv[]
) {
    double rho = 0;
    int points = 10;
    double increment = 1.0/points;
    opdXcosm<double> opdX(
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
	160,
	160,
        0.9,            // numerical aperture
        0.000488        // fluorescence wavelenght
    );
    opdGibsonLani<double> opdG(
        0.050,          // specimen thickness 
        0.16,           // immersion thickness design
        0.16,           // immersion thickness actual
        0.170,          // coverglass thickness design
        0.120,          // coverglass thickness actual
        1.33,           // specimen refractive index
        1.0,            // immersion refractive index design
        1.0,            // immersion refractive index actual
        1.525,          // coverglass refractive index design
        1.525           // coverglass refractive index actual
    );
    opdTorokVarga<double> opdT(
        0.050,          // specimen thickness 
        0.16,           // immersion thickness design
        0.16,           // immersion thickness actual
        0.170,          // coverglass thickness design
        0.120,          // coverglass thickness actual
        1.33,           // specimen refractive index
        1.0,            // immersion refractive index design
        1.0,            // immersion refractive index actual
        1.525,          // coverglass refractive index design
        1.525           // coverglass refractive index actual
    );
    for ( int i = 0; i < points; i++ ) 
    {
	rho += increment;
	std::cout << rho << ", "<< opdX(rho) <<","<< opdG(rho)<<", "<<opdT(rho)<< std::endl;	
    }
    return 0;
}
