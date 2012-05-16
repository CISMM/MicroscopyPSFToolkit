
/* Copyright 1993 Washington University Biomedical Computer Laboratory */

#include <stdio.h>
#include <math.h>

#define PI 3.1415926535979

#define ReDat(i) *(rdata + i)
#define ImDat(i) *(idata + i)
#define SWAP(a,b) { tempx = a; a = b; b = tempx; }

void fftn(float *rdata, float *idata, unsigned nn[], int ndim, int isign)
{
int idim;
unsigned i1, i2rev, i3rev, ibit;
unsigned ip2, ifp1, ifp2, k2, n;
unsigned nprev = 1, nrem, ntot = 1;
register unsigned i2, i3;
float *ptrr,*ptri;
double theta;
float wx,wy,wpx,wpy,fntot;
double wtemp;
float tempx,tempy,wtx,wty;
double t1, t2;

/*      Compute total number of complex values  */
for (idim = 0; idim < ndim; ++idim)
	ntot *= nn[idim];
fntot = (float)ntot;
for (idim = ndim - 1; idim >= 0; --idim) {
	n = nn[idim];

	nrem = ntot / (n * nprev);
	ip2 = nprev * n;        /*      Unit step for next dimension */
	i2rev = 0;              /*      Bit reversed i2 */

	/*      This is the bit reversal section of the routine */
	/*      Loop over current dimension     */
	for (i2 = 0; i2 < ip2; i2 += nprev) {
		if (i2 < i2rev)
		/*      Loop over lower dimensions      */
			for (i1 = i2; i1 < i2 + nprev; ++i1)
			/*      Loop over higher dimensions  */
 			   for (i3 = i1; i3 < ntot; i3 += ip2) {
				i3rev = i3 + i2rev - i2;
				SWAP(ReDat(i3),ReDat(i3rev));
				SWAP(ImDat(i3),ImDat(i3rev));
				}

			ibit = ip2;
			/*      Increment from high end of i2rev to low */
			do {
				ibit >>= 1;
				i2rev ^= ibit;
			} while (ibit >= nprev && !(ibit & i2rev));
		}

/*      Here begins the Danielson-Lanczos section of the routine */
		/*      Loop over step sizes    */
		for (ifp1 = nprev; ifp1 < ip2; ifp1 <<= 1) {
			ifp2 = ifp1 << 1;
			/*  Initialize for the trig. recurrence */
			theta = isign * 2.0 * PI / (ifp2 / nprev);
			wpx = sin(0.5 * theta);
			wpx *= -2.0 * wpx;
			wpy = sin(theta);
			wx = 1.0;
			wy = 0.0;

		/*  Loop by unit step in current dimension  */
		for (i3 = 0; i3 < ifp1; i3 += nprev) {
			/*      Loop over lower dimensions      */
			for (i1 = i3; i1 < i3 + nprev; ++i1)
				/*  Loop over higher dimensions */
				for (i2 = i1; i2 < ntot; i2 += ifp2) {
					/*      Danielson-Lanczos formula */
					k2 = i2 + ifp1;
					wtx = ReDat(k2);
					wty = ImDat(k2);

/*	Complex multiply using 3 real multiplies.  Should usually be faster.	*/
ReDat(k2) = ReDat(i2) - (tempx = (t1=wx*wtx) - (t2=wy*wty));
ImDat(k2) = ImDat(i2) - (tempy = (wx + wy) * (wtx + wty) - t1 - t2);
					ReDat(i2) += tempx;
					ImDat(i2) += tempy;
					}
				/*      Trigonometric recurrence        */
				wtemp = wx;
/*	Complex multiply using 3 real multiplies.	*/
				wx += (t1 = wx * wpx) - (t2 = wy * wpy);
				wy += (wtemp + wy) * (wpx + wpy) - t1 - t2;
			}
		}
	nprev *= n;
	}
/* Normalize if reverse */
if(isign == 1)
{
  ptrr = rdata; 
  ptri = idata;
  for(idim=0;idim<ntot;idim++)
  {
     *ptrr /= fntot;
     ptrr++;
     *ptri /= fntot;
     ptri++;
  }
}
}
