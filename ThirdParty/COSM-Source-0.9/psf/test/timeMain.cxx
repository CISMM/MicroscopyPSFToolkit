#include "psf/qtrapIntegrator.h"
#include "psf/qsimpIntegrator.h"
#include "psf/qrombIntegrator.h"
#include "psf/radialPSF.h"
#include "psf/gibsonLaniFunctor.h"
#include "psf/opdXcosm.h"

#include <blitz/timer.h>
#include <iostream> 
#include <stdio.h>

using namespace cosm;

#define REAL double

int main ( int argc, char* argv[] ) {

    REAL na = 0.9;
    REAL lambda = 488*1e-6;
    REAL k = M_PI*2.0/lambda;
    REAL ns = 1.33;
    REAL eps = 1e-4;
    int jmax = 20;

    opdXcosm<REAL> opd(
        0.050,          // specimen thickness
        0.16,           // immersion thickness design
        0.16,           // immersion thickness actual
        0.170,          // coverglass thickness design
        0.120,          // coverglass thickness actual
        ns,             // specimen refractive index
        1.0,            // immersion refractive index design
        1.0,            // immersion refractive index actual
        1.525,          // coverglass refractive index design
        1.525,          // coverglass refractive index actual
        160,            // tube length design
        160,            // tube length actual
        100.0,          // lateral magnification
        na              // numerical aperture (NA)
    );
    const int max_int = 3;
    GibsonLaniFunctor<REAL> gibsonLani(k, opd);
    QTrapIntegrator<REAL> qtrap(gibsonLani, 0,ns*ns, eps, jmax);
    QSimpIntegrator<REAL> qsimp(gibsonLani, 0,ns*ns, eps, jmax);
    QRombIntegrator<REAL> qromb(gibsonLani, 0,ns*ns, eps, jmax);
    Integrator<REAL>* integrator[max_int] = {&qtrap, &qsimp, &qromb};
    RadialPSF<REAL> radialPSF(
        8,              // size of the image in x and y
        8,              // size of the image in z
        0.068*1e-3,     // pixel size in z and y of the image (mm)
        0.1*1e-3,       // pixel size in z of the image (mm)
        &gibsonLani,
	integrator[0],
        lambda,         // fluorescence wavelenght (lambda) (mm)
        na,             // numerical aperture (NA)
        RadialPSF<REAL>::OPTICAL_SECTIONING_WIDEFIELD
    );

    blitz::Timer timer;
    for ( int i = 2; i < max_int; i++ ) 
    {
	radialPSF.setIntegrator( integrator[i] );
	timer.start();
        radialPSF.evaluate();
	timer.stop();
        printf("time for integration: %e\n", timer.elapsedSeconds());
    }
    return 0;
}
