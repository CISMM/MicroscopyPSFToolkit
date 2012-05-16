/***************************************************************************
  COPYRIGHT 1996 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
/* 

REAL_FFT3D.C

	FFT routines for 3-Dimensional real data.  If working on a
Silicon Graphics with the libcomplib.sgimath.a library, use the
routines developed especially for the SG cache and hardware (3X faster), 
else use simple 3-pass DFT. Code for 1-D DFT (fftl) and 1-D REAL DFT (fftr)
are from Numerical Recipes in C.

   NOTE: Input data must be stored in REAL array of size (NX+2,NY,NZ)
	 with the original data stored in the first NX columns. Last two
	 columns are for the sfft3du/fftr routine.
 
	 Input array is over-written with output. For REVERSE, last two
	 columns are meaningless.

Author:
	Keith Doolittle
	Routines fftr() and fftl() from Numerical Recipes in C
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft3d.h"

#define MAX_DIM 16384		/* Maximum size of one dimension of data */

#ifdef SGI_FFT

int   	done_init = 0;			/* Initialized yet?   	*/
int   	init_nx,init_ny,init_nz;	/* Dims done by init 	*/
float 	*work;				/* Initialization array */

int	done_init1D = 0;
int	init_nx1D;
fcomplex	work1D[MAX_DIM+15];		/* 1D Initialization array */

void real_fft3d(float *x,int nx, int ny, int nz, int xdim, int ydim, int dir)
{
if((done_init == 0)||(nx != init_nx)||(ny != init_ny)||(nz != init_nz)) {
 if(done_init != 0) free((char*)work);
 if((work=(float*)malloc(sizeof(float)*(100+45+nx+2*ny+2*nz)))==(float*)NULL) {
    fprintf(stderr,"ERROR: Out of memory in real_fft3d (SGI)\n");
    exit(1);
    }
  work = sfft3dui(nx,ny,nz,work);
  init_nx = nx; init_ny = ny; init_nz = nz; done_init = 1;
  }
(void)sfft3du(dir,nx,ny,nz,x,xdim,ydim,work);
if(dir == REVERSE)
  (void)sscal3d(nx,ny,nz,1.0/((float)(nx*ny*nz)),x,xdim,ydim);
}

/* complex -- complex 1-D FFT */
void cplx_fft1d(fcomplex *x,int nx, int dir)
{
if((done_init1D == 0)||(nx != init_nx1D)) {
   init_nx1D = nx;
   (void)cfft1di(nx,work1D);
   }
(void)cfft1d(dir,nx,x,1,work1D);
}


#else /* ifdef SGI_FFT */

#define Dat(dat,x,y,z,nx,nxy) *(dat + (x) + (y)*(nx) + (z)*(nxy))
#define Ptr(dat,x,y,z,nx,nxy)  (dat + (x) + (y)*(nx) + (z)*(nxy))


fcomplex wk[MAX_DIM+1];		/* Working array 			*/

void fftr(int,int);		/* Real    --> Complex FFT 1-D          */
void fftl(int,int);             /* Complex --> Complex FFT 1-D          */


/* Complex conjugate */
fcomplex conjg(fcomplex x)
{
fcomplex fc;
fc.re = x.re;
fc.im = -x.im;
return(fc);
}
  
/* complex-->complex 1-D FFT */
void cplx_fft1d(fcomplex *x, int nx, int dir)
{
int i;
fcomplex *f = x;

for(i=0;i<nx;i++,f++) 
  wk[i] = *f;
fftl(nx,dir);
}


void scalefft(float *xx,int nx, int ny, int nz, int dx, int dy)
{
int x,y,z,deltx,pl_size;
float fscl,*fptr;

pl_size = dx*dy;
deltx   = dx-nx;
fscl    = 1.0/((float)(nx*ny*nz/2));
for(z=0;z<nz;z++) {
  fptr = xx + z*pl_size;
  for(y=0;y<ny;y++,fptr += deltx)
    for(x=0;x<nx;x++,fptr++)
        *fptr = *fptr * fscl;
  }
}

/***********************************************************************
   Three pass 3-D DFT. Order = XYZ (FORWARD) or ZYX (REVERSE). For
   X dimension, use fftr which takes input array of size N+2 and
   output's complex array of size (N+2)/2 (for REAL sequences). Others
   use standard fftl complex DFT.
***********************************************************************/

void real_fft3d(float *xx,int nx, int ny, int nz, int xdim, int ydim, int dir)
{
int cx,cxy,Halfnx;
int i1,i2,i3,xadd,pl_size;
float *f,fscl;
fcomplex *cptr,fcomp;
fcomplex *x = (fcomplex*)xx;

if((nx <0)||(nx > MAX_DIM)||(ny <0)||(ny > MAX_DIM)||(nz <0)||(nz > MAX_DIM)){
  fprintf(stderr,"ERROR: real_fft3d recompile for larger than %d ",MAX_DIM);
  fprintf(stderr,"in one dimension (DIMS %d,%d,%d)\n",nx,ny,nz);
  exit(1);
  }

/* old -- before xdim/ydim add
Halfnx = nx/2;
cx = Halfnx+1;
cxy = cx*ny;
*/
Halfnx = nx/2;
cx     = xdim/2;
cxy    = cx*ydim;
xadd   = xdim-nx;
pl_size= xdim*ydim;

if(dir == FORWARD) {

for(i3=0;i3<nz;i3++)
 for(i2=0;i2<ny;i2++) {
  for(i1=0;i1<Halfnx;i1++)
	wk[i1] = Dat(x,i1,i2,i3,cx,cxy);
  fftr(Halfnx,dir);
  cptr = Ptr(x,0,i2,i3,cx,cxy);
  cptr->re = wk[0].re; cptr->im = 0.0;
  cptr = Ptr(x,Halfnx,i2,i3,cx,cxy);
  cptr->re = wk[0].im; cptr->im = 0.0;
  for(i1=1;i1<Halfnx;i1++)
        Dat(x,i1,i2,i3,cx,cxy) = wk[i1];
  }

for(i3=0;i3<nz;i3++)
 for(i1=0;i1<=Halfnx;i1++) {
  for(i2=0;i2<ny;i2++)
	wk[i2] = Dat(x,i1,i2,i3,cx,cxy);
  fftl(ny,dir);
  for(i2=0;i2<ny;i2++)
	Dat(x,i1,i2,i3,cx,cxy) = wk[i2];
  }

if(nz == 1) return;

for(i1=0;i1<=Halfnx;i1++)
 for(i2=0;i2<ny;i2++) {
  for(i3=0;i3<nz;i3++)
	wk[i3] = Dat(x,i1,i2,i3,cx,cxy);
  fftl(nz,dir);
  for(i3=0;i3<nz;i3++) 
	Dat(x,i1,i2,i3,cx,cxy) = wk[i3];
  }


}
else { /* dir == REVERSE */

if(nz == 1) goto SkipZ;
for(i1=0;i1<=Halfnx;i1++) 
 for(i2=0;i2<ny;i2++) {
  for(i3=0;i3<nz;i3++)
	wk[i3] = Dat(x,i1,i2,i3,cx,cxy);
  fftl(nz,dir);
  for(i3=0;i3<nz;i3++)
	Dat(x,i1,i2,i3,cx,cxy) = wk[i3];
  }

SkipZ:
for(i3=0;i3<nz;i3++)
 for(i1=0;i1<=Halfnx;i1++) {
  for(i2=0;i2<ny;i2++)
	wk[i2] = Dat(x,i1,i2,i3,cx,cxy);
  fftl(ny,dir);
  for(i2=0;i2<ny;i2++)
	Dat(x,i1,i2,i3,cx,cxy) = wk[i2];
  }

for(i3=0;i3<nz;i3++)
 for(i2=0;i2<ny;i2++) {
  cptr = Ptr(x,0,i2,i3,cx,cxy);
  wk[0].re = cptr->re;
  cptr = Ptr(x,Halfnx,i2,i3,cx,cxy);
  wk[0].im = cptr->re;
  for(i1=1;i1<=Halfnx;i1++)
	wk[i1] = Dat(x,i1,i2,i3,cx,cxy);
  fftr(Halfnx,dir);
  for(i1=0;i1<Halfnx;i1++)
	Dat(x,i1,i2,i3,cx,cxy) = wk[i1];
  }
scalefft(xx,nx,ny,nz,xdim,ydim);

 }
}


/* Takes FFT of real data which is 2N length */

void fftr(int n, int direc)
{
float *dat,*fptr;
int i,i1,i2,i3,i4,n2p3;
float c1,c2,h1r,h1i,h2r,h2i;
double wr,wi,wpr,wpi,wtemp,theta;

fptr = (float*)&wk[0];
dat = fptr-1;

c1 = 0.5;
theta = 3.141592653589793/(double)n;
if(direc == FORWARD) {
  c2 = -0.5; fftl(n,1);
  }
else {
  c2 = 0.5; theta = -theta;
  }

wtemp = sin(0.5*theta);
wpr = -2.0*wtemp*wtemp;
wpi = sin(theta);
wr = 1.0+wpr; wi = wpi;
n2p3 = 2*n+3;
for(i=2;i<=n/2;i++) {
  i4 = 1+(i3=n2p3-(i2=1+(i1=i+i-1)));
  h1r = c1*(dat[i1]+dat[i3]);
  h1i = c1*(dat[i2]-dat[i4]);
  h2r = -c2*(dat[i2]+dat[i4]);
  h2i = c2*(dat[i1]-dat[i3]);
  dat[i1] = h1r+wr*h2r-wi*h2i;
  dat[i2] = h1i+wr*h2i+wi*h2r;
  dat[i3] = h1r-wr*h2r+wi*h2i;
  dat[i4] = -h1i+wr*h2i+wi*h2r;
  wr = (wtemp=wr)*wpr-wi*wpi+wr;
  wi = wi*wpr+wtemp*wpi+wi;
  }
if(direc == FORWARD) {
  dat[1] = (h1r=dat[1])+dat[2];
  dat[2] = h1r-dat[2];
  }
else {
  dat[1] = c1*((h1r=dat[1])+dat[2]);
  dat[2] = c1*(h1r-dat[2]);
  fftl(n,REVERSE);
  }
}

#define SWAP(a,b) tempr = (a); (a)=(b); (b)=tempr;

void fftl(int nn, int direc)
{
float *dat,*fptr;
int n,mmax,m,j,istep,i;
double wtemp,wr,wpr,wpi,wi,theta;
float tempr,tempi;

/* This routine is 1 based, so start at wk[0]-sizeof(float) */
fptr = (float*)&wk[0];
dat = fptr-1;
n = nn<<1;
j = 1;
for(i=1;i<n;i+=2) {
 if(j>i) {
  SWAP(dat[j],dat[i]);
  SWAP(dat[j+1],dat[i+1]);
  }
  m = n>>1;
  while(m >= 2 && j > m) {
   j -= m; m >>=1;
   }
  j += m;
}
mmax = 2;
while (n > mmax) {
  istep = 2*mmax;
  theta = 6.28318530717959/(direc*mmax);
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  wr = 1.0;
  wi = 0.0;
  for(m=1;m<mmax;m+=2) {
   for(i=m;i<=n;i+=istep) {
     j=i+mmax;
     tempr = wr*dat[j]-wi*dat[j+1];
     tempi = wr*dat[j+1]+wi*dat[j];
     dat[j] = dat[i]-tempr;
     dat[j+1] = dat[i+1]-tempi;
     dat[i] += tempr; dat[i+1] += tempi;
   }
   wr = (wtemp=wr)*wpr-wi*wpi+wr;
   wi = wi*wpr+wtemp*wpi+wi;
  }
  mmax=istep;
}
}

#endif /* ifdef SGI_FFT */
