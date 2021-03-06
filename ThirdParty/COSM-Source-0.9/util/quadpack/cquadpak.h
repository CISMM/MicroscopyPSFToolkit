/* [einirv] - changed "double f()" to "double f(double)" 
 *            added __cplusplus idfdefs
 */
#ifndef _CQADPAK_H
#define _CQADPAK_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* Integration routines */
/* Gauss-Kronrod for integration over finite range. */
double G_K15(double f(double, void *),void* cbData,double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K21(double f(double, void *),void* cbData,double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K31(double f(double, void *),void * cbData,double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K41(double f(double, void *),void * cbData,double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K51(double f(double, void *),void * cbData,double a,double b,double *abserr,
	double *resabs, double *resasc);
double G_K61(double f(double, void *),void * cbData,double a,double b,double *abserr,
	double *resabs, double *resasc);

/* Gauss-Kronrod for integration over infinite range. */
double G_K15I(double f(double, void *),void * cbData,double boun,int inf,double a,double b,
	double *abserr,double *resabs, double *resasc);

/* Gauss-Kronrod for integration of weighted function. */
double G_K15W(double f(double, void *),void * cbData,double w(),double p1,double p2,double p3,
	double p4,int kp,double a,double b,double *abserr,
	double *resabs, double *resasc);
double dqext(int *n,double epstab[],double *abserr,
	double res3la[],int *nres);
void dqsort(int limit,int last,int *maxerr,double *ermax,
	double elist[],int iord[],int *nrmax);
double dqagi(double f(double, void*),void * cbData,double bound,int inf,double epsabs,
	double epsrel,double *abserr,int *neval,int *ier);
double dqags(double f(double, void*),void * cbData,double a,double b,double epsabs,
	double epsrel,double *abserr,int *neval,int *ier);
double dqagp(double f(double, void *),void * cbData,double a,double b,int npts2,double *points,
	double epsabs,double epsrel,double *abserr,int *neval,int *ier);
double dqng(double f(double, void *),void * cbData,double a,double b,double epsabs,double epsrel,
	double *abserr,int *neval,int *ier);
double dqag(double f(double, void *),void *cbData,double a,double b,double epsabs,double epsrel,
	int irule,double *abserr,int *neval,int *ier);
double dqage(double f(double, void *),void *cbData,double a,double b,double epsabs,double epsrel,
        int irule,double *abserr,int *neval,int *ier,int *last);
double dqwgtc(double x,double c,double p2,double p3,double p4,
	int kp);
double dqwgto(double x,double omega,double p2,double p3,double p4,
	int integr);
double dqwgts(double x,double a,double b,double alpha,double beta,
	int integr);
void dqcheb(double *x,double *fval,double *cheb12,double *cheb24);
double dqc25o(double f(double, void *),void * cbData,double a,double b,double omega,int integr,
	int nrmom,int maxp1,int ksave,double *abserr,int *neval,
	double *resabs,double *resasc,int *momcom,double **chebmo);		
double dqfour(double f(double, void *),void * cbData,double a,double b,double omega,int integr,
    double epsabs,double epsrel,int icall,int maxp1,
    double *abserr,int *neval,int *ier,int *momcom,
	double **chebmo);
double dqawfe(double f(double, void *),void * cbData,double a,double omega,int integr,double epsabs,
    int limlst,int maxp1,double *abserr,int *neval,int *ier,
    double *rslst,double *erlist,int *ierlst,double **chebmo);
double dqawf(double f(double, void *),void * cbData,double a,double omega,int integr,double epsabs,
	double *abserr,int *neval,int *ier);
double dqawo(double f(double, void *),void * cbData,double a,double b,double omega,int integr,double epsabs,
	double epsrel,double *abserr,int *neval,int *ier);
double dqaws(double f(double, void *),void * cbData,double a,double b,double alfa,double beta,int wgtfunc,
    double epsabs,double epsrel,double *abserr,int *neval,int *ier);
double dqawse(double f(double, void *),void * cbData,double a,double b,double alfa,double beta,
    int wgtfunc,double epsabs,double epsrel,double *abserr,
    int *neval,int *ier);
void dqmomo(double alfa,double beta,double ri[],double rj[],double rg[],
    double rh[],int wgtfunc);
double dqc25s(double f(double, void *),void * cbData,double a,double b,double bl,double br,double alfa,
    double beta,double ri[],double rj[],double rg[],double rh[],
    double *abserr,double *resasc,int wgtfunc,int *nev);
double dqc25c(double f(double, void *),void * cbData,double a,double b,double c,double *abserr,
    int *krul,int *neval);
double dqawc(double f(double, void *),void * cbData,double a,double b,double c,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier);
double dqawce(double f(double, void *),void * cbData, double a,double b,double c,double epsabs,
    double epsrel,double *abserr,int *neval,int *ier);

double G_B15(double f(double, void *),void * cbData,double a,double b,double *abserr,
	double *resabs, double *resasc);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* _CQUADPAK_H */
