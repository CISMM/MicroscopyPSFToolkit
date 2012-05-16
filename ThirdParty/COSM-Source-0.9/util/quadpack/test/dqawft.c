#include <stdio.h>
#include <float.h>
#include "cquadpak.h"
#include "dqdefs.h"
#include <math.h>

double f(double x)
{
    return ((x > 0.0) ? (1.0/sqrt(x)) : 0.0);
}
double fx(double x)
{
	return exp(-x);
}

/* [einirv] - change main return type from void to int */
int main()
{
	double a,omega,result,abserr,epsabs;
	int ier,neval;

	a = 0.0;
    omega = 0.5 * Pi;
    epsabs = 1.0e-8;
    result =dqawf(fx,a,omega,COSINE,epsabs,&abserr,
		&neval,&ier);
	printf("\nresult = %.17lg\n",result);
	printf("abserr = %.17lg\n",abserr);
	printf("neval = %d\n",neval);
	return 0;
}
