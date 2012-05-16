/*************************************************************
    spline.c -- translated by f2c (version 19980913).
               obtained from netlib (URL: www.netlib.org)
    
    The array y() of samples to be used for the interpolation 
    was changed from double to float by C. Preza
*************************************************************/
/*  METHOD */
/*  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed */
/*  for a cubic interpolating spline */

/*    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3 */

/*    for  x(i) .le. x .le. x(i+1) */

/*  input.. */

/*    n = the number of data points or knots (n.ge.2) */
/*    x = the abscissas of the knots in strictly increasing order */
/*    y = the ordinates of the knots */

/*  output.. */

/*    b, c, d  = arrays of spline coefficients as defined above. */

/*  using  p  to denote differentiation, */

/*    y(i) = s(x(i)) */
/*    b(i) = sp(x(i)) */
/*    c(i) = spp(x(i))/2 */
/*    d(i) = sppp(x(i))/6  (derivative from the right) */
/*  the accompanying function subprogram  seval  can be used */
/*  to evaluate the spline. */

#include "f2c.h"

/* Subroutine */ 
int spline_(n, x, y, b, c__, d__)
integer *n; 
doublereal *x; 
float *y; 
doublereal *b, *c__, *d__;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal t;
    static integer ib, nm1;


    /* Parameter adjustments */
    /* decrement the array pointers by one because the loop index starts
       at 1 instead of 0 */
    --d__;
    --c__;
    --b;
    --y;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    if (*n < 2) {
	return 0;
    }
    if (*n < 3) {
	goto L50;
    }

/*  set up tridiagonal system */

/*  b = diagonal, d = offdiagonal, c = right hand side. */

    d__[1] = x[2] - x[1];
    c__[2] = (y[2] - y[1]) / d__[1];
    i__1 = nm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	d__[i__] = x[i__ + 1] - x[i__];
	b[i__] = (d__[i__ - 1] + d__[i__]) * (float)2.;
	c__[i__ + 1] = (y[i__ + 1] - y[i__]) / d__[i__];
	c__[i__] = c__[i__ + 1] - c__[i__];
/* L10: */
    }

/*  end conditions.  third derivatives at  x(1)  and  x(n) */
/*  obtained from divided differences */

    b[1] = -d__[1];
    b[*n] = -d__[*n - 1];
    c__[1] = (float)0.;
    c__[*n] = (float)0.;
    if (*n == 3) {
	goto L15;
    }
    c__[1] = c__[3] / (x[4] - x[2]) - c__[2] / (x[3] - x[1]);
    c__[*n] = c__[*n - 1] / (x[*n] - x[*n - 2]) - c__[*n - 2] / (x[*n - 1] - 
	    x[*n - 3]);
/* Computing 2nd power */
    d__1 = d__[1];
    c__[1] = c__[1] * (d__1 * d__1) / (x[4] - x[1]);
/* Computing 2nd power */
    d__1 = d__[*n - 1];
    c__[*n] = -c__[*n] * (d__1 * d__1) / (x[*n] - x[*n - 3]);

/*  forward elimination */

L15:
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	t = d__[i__ - 1] / b[i__ - 1];
	b[i__] -= t * d__[i__ - 1];
	c__[i__] -= t * c__[i__ - 1];
/* L20: */
    }

/*  back substitution */

    c__[*n] /= b[*n];
    i__1 = nm1;
    for (ib = 1; ib <= i__1; ++ib) {
	i__ = *n - ib;
	c__[i__] = (c__[i__] - d__[i__] * c__[i__ + 1]) / b[i__];
/* L30: */
    }

/*  c(i) is now the sigma(i) of the text */

/*  compute polynomial coefficients */

    b[*n] = (y[*n] - y[nm1]) / d__[nm1] + d__[nm1] * (c__[nm1] + c__[*n] * (
	    float)2.);
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = (y[i__ + 1] - y[i__]) / d__[i__] - d__[i__] * (c__[i__ + 1] 
		+ c__[i__] * (float)2.);
	d__[i__] = (c__[i__ + 1] - c__[i__]) / d__[i__];
	c__[i__] *= (float)3.;
/* L40: */
    }
    c__[*n] *= (float)3.;
    d__[*n] = d__[*n - 1];
    return 0;

L50:
    b[1] = (y[2] - y[1]) / (x[2] - x[1]);
    c__[1] = (float)0.;
    d__[1] = (float)0.;
    b[2] = b[1];
    c__[2] = (float)0.;
    d__[2] = (float)0.;
    return 0;
} /* spline_ */

