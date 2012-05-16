/***************************************************************************
 COPYRIGHT 2000 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#ifndef WIN32
#include <sys/param.h>
#include <sys/time.h>
#endif
#include "washu.h"

extern float deltaxy;
extern float deltaxy_nyq;
extern char lognm[];
extern void WritePlane(int plane, float *outvol, float *inslice,
                        int w, int h, int inx ,int iny);
extern int spline_(int *points, double *coord, float *samples, double *b_coef,
double *c_coef, double *d_coef);

#define DatOf1D(dat,x)                  *(dat + (x))
#define DatOf2D(dat,x,y,wid)            *(dat + (x) + (y)*wid)
#define MAXPATHLEN 256

#ifdef DEBUG
extern void ShowValues(float *d, int w, int h, char *fname, int tpe);

#endif

void Save2D(float *dat, int N, char *fname)
{
osm_ds head;
memset((void*)&head,0,sizeof(osm_ds));
head.data = (char*)dat;
head.mode = WU_FLOAT;
head.user15 = TVAL;
head.user16_footer_size = 0;
head.nx = N;
head.ny = N;
head.nz = 1;
write_dataset(&head,fname);
}

/* sums at least 5 neighboring pixels to obtain a more accurate pixel for
an undersampled PSF*/

void Sum4N4NToNN(float *src, float *dest, int N, int osamp)
{
int i,j,jj,ii,jind,iind;

int N4 = osamp*N;
int mo2 = osamp/2;
double dsum;


for(j=0;j<N;j++)  
  for(i=0;i<N;i++) {

	dsum = 0.0;
	for(jj=osamp*j-mo2;jj<=osamp*j+mo2;jj++) {
	  if(jj < 0) jind = -jj;
	  else if(jj >= N4) jind = 2*N4 - jj - 2;
	  else jind = jj;
	  for(ii=osamp*i-mo2;ii<=osamp*i+mo2;ii++) {
	  	if(ii < 0) iind = -ii;
	  	else if(ii >= N4) iind = 2*N4 - ii - 2;
		else iind = ii;

		dsum += DatOf2D(src,iind,jind,N4);
		}
	  }

	DatOf2D(dest,i,j,N) = dsum;
	}
}

/*************************************************************************************
  ROTSUM :  PSF volume creation for widefield or 2-photon undersampled PSFs.
	    This routine reads in the radial section of the PSF that was
	    generated by xcosm_psf and radially sweeps it 360 degrees
	    (centered at 0,0) by interpolating the values.
	    Because the input pixel size is larger than the Nyquist rate,
	    the values are interpolated at a finer resolution than the
	    input pixel size in order to avoid problems with undersampling.
	    The final PSF values are computed by summing at least 5
	    neighboring pixel values.
	    The program also checks for a symmetric psf, in which case
	    the above procedure need only be performed for N/2+1 rows	     
	    and then duplicated.
	    The resulting psf is centered at 0,0,0.	

 Author: Keith Doolittle
*************************************************************************************/

void rotsum(float *outimg, osm_ds *head)
{
char	time_str[MAXPATHLEN];
char	temp_str[MAXPATHLEN];
FILE	*fplog;
int   	oversmp,symmetric;
int   	N,N2,N4,Nxy,Nz,Nr,Nzwanted;
float 	deltaxy1,deltar;
float 	*psf4N4N,*psfNN,*outpsf,*rowdat;
int   	ix,iy,iz,ir;
int     osamp;
float 	r_hat,r,x,y;
float 	interp_val;
double  cord, sq_cord;
double  *bcf, *ccf, *dcf, *dr_cord;

time_t 	clock;

outpsf    = (float*)head->data;		/* XY slice data */
Nz 	  = head->ny;
Nr	  = head->nx;

oversmp   = head->xstart;		/* oversmp */
Nxy       = head->ystart;
symmetric = head->zstart;		/* symmetr */
deltar    = head->xlength;
deltaxy   = head->ylength;



/* oversampling rate    */
/* make it beat nyquist */
osamp = (int)(deltaxy/deltaxy_nyq) + 1;

if(osamp < 5) osamp = 5;

/* make odd */
if(!(osamp & 1)) osamp++;

/* New deltaxy for oversampling */
deltaxy1 = deltaxy/((float)osamp);

#ifdef DEBUG
fprintf(stderr,"rotsum; DXY: %.4f DXY_NYQ: %.4f SAMP = %d NEW_DXY: %.4f\n",
		deltaxy,deltaxy_nyq,osamp,deltaxy1);
#endif

N      = Nxy;			
N4     = N*osamp; 
N2     = N4/2; 

if(symmetric == 0)
   Nzwanted = Nz;
else
   Nzwanted = Nz/2+1;

/* Complex psf (4N+2)/2 x 4N oversampled 2-D plane */

if((psf4N4N=(float*)calloc(sizeof(float),N4*N4))==(float*)NULL) {
 fprintf(stderr,"rotsum; ERROR: Out of memory (2).\n");
 exit(1);
 }

/* Complex psf slice (N+2)/2 x N 2-D plane */

if((psfNN=(float*)calloc(sizeof(float),N*N))==(float*)NULL) {
 fprintf(stderr,"rotsum; ERROR: Out of memory (6).\n");
 exit(1);
 }

rowdat = (float*)outpsf;

if((dr_cord=(double*)calloc(sizeof(double),Nr))==(double*)NULL) {
 fprintf(stderr,"rotsum; ERROR: Out of memory (6).\n");
 exit(1);
}
if((bcf=(double*)calloc(sizeof(double),Nr))==(double*)NULL) {
 fprintf(stderr,"rotsum; ERROR: Out of memory (6).\n");
 exit(1);
}
if((ccf=(double*)calloc(sizeof(double),Nr))==(double*)NULL) {
 fprintf(stderr,"rotsum; ERROR: Out of memory (6).\n");
 exit(1);
}
if((dcf=(double*)calloc(sizeof(double),Nr))==(double*)NULL) {
 fprintf(stderr,"rotsum; ERROR: Out of memory (6).\n");
 exit(1);
}

/* initialize radial coordinates array for the spline call */
 for(ir=0;ir<Nr;ir++)
	 dr_cord[ir] = ir * deltar;

for(iz=0;iz<Nzwanted;iz++) {

/* compute the cubic spline interpolation coefficients: bcf, ccf, dcf */
spline_(&Nr,dr_cord, rowdat,bcf,ccf,dcf);

 /* Sweep psf radial data into 1 plane */
 for(iy=0;iy<=N2;iy++) {
  y = (float)iy * deltaxy1;
  for(ix=0;ix<=N2;ix++) {
   x = (float)ix * deltaxy1;
   r = (float)sqrt(x*x + y*y);
   r_hat = r/deltar;

   ir = (int)r_hat;
   cord = r - dr_cord[ir];
   sq_cord = cord * cord;

/* interpolated value based on cubic spline interpolation equation */
   interp_val = rowdat[ir] + bcf[ir]*cord + ccf[ir]*sq_cord + dcf[ir]*cord*sq_cord;

   /* Upper left */
   if((iy < N2)&&(ix < N2))
     DatOf2D(psf4N4N,ix,iy,N4) = interp_val;
   /* Upper right */
   if(ix > 0)
     DatOf2D(psf4N4N,N4-ix,iy,N4) = interp_val;
   /* Lower left */
   if(iy > 0)
     DatOf2D(psf4N4N,ix,N4-iy,N4) = interp_val;
   /* Lower right */
   if((ix > 0)&&(iy > 0))
     DatOf2D(psf4N4N,N4-ix,N4-iy,N4) = interp_val;
   }
  }

#ifdef DEBUG
if(iz == 0) {
   //ShowValues(psf4N4N,N4,N4,"orig4N",4);
   Save2D(psf4N4N,N4,"psf_slice4N");
   }
#endif

 /*------------------------------------------------------------------
   Now, sum in X/Y direction by at least 5 pixels in each direction
   ------------------------------------------------------------------*/

 Sum4N4NToNN(psf4N4N,psfNN,N,osamp);

#ifdef DEBUG
if(iz == 0) {
   //ShowValues(psfNN,N,N,"origN",4);
   Save2D(psfNN,N,"psf_sliceN");
   }
#endif

 /* Copy to 3-D PSF */

 WritePlane(iz,outimg,psfNN,N,N,N,N);
 if((symmetric != 0)&&(iz > 0)&&(iz < Nzwanted-1))
	WritePlane(Nz-iz,outimg,psfNN,N,N,N,N);

 /* Increment row counter to point to next row */
 rowdat += Nr;

 if((fplog=fopen(lognm,"a"))==(FILE*)NULL)
   fprintf(stderr,"rotsum; WARNING: Can't open logfile `%s' for write/append (rotnormal)\n",lognm);
 else {
      clock = time(NULL);
      strcpy(time_str,asctime(localtime(&clock)));
      time_str[strlen(time_str)-1] = '\0';
      sprintf(temp_str," Completed plane %d",iz);
      strcat(time_str,temp_str);
      fprintf(fplog,"rotsum; %s\n",time_str);
      fclose(fplog);
      }
 }

}

