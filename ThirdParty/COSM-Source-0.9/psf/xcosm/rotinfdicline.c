/*
	INF_DIC_PSF.C       : program that computes the coherent 2D psf for DIC
	Last modified by: Chrysanthe Preza on 11/96
	                  to accomodate the DIC optics
	    Washington University Biomedical Computer Laboratory.

	this program is to be used with the rotation so note that
	the image is 2 times larger than the desired one 
	but because we want the same sampling, a factor of 2 is used
	in determining the spacing in frequency.
*/
#include <stdio.h>
#ifndef WIN32
#include <sys/param.h>
#endif
#include <math.h>
#include "washu.h"

#define MAXPATHLEN 256
#define PI      3.1415926535
#define DatOf1D(dat,x)			*(dat + (x))
#define DatOf2D(dat,x,y,wid) 		*(dat + (x) + (y)*wid)

/* use 0.003400 mm for the 10x0.3 NA lens */
#define sampling_xy	0.003400  /*[mm] in image space   */

/* use 5.3e-5  [mm] for the case used by Tim Holmes */
/*#define sampling_xy	5.3e-5  [mm] in image space; this is for T. Holmes
				set up. The mag=1.0 so dxy is same in object
				and image space */

int N, Nz;

void rotinfdicline(float *outpsf_re,
                float * outpsf_im,
                osm_ds *head_re,
				osm_ds *head_im,
                int bin,
                float DeltaX,
                float DeltaPhi,
                float AmpRatio,
                int tandem)
{
/* values for DIC*/
int num_pix;
float twoPInum_pix;
float sinfactor, cosfactor, ratfactor, angle;
float f, deltaf, delta, fmax, fmax_half;

FILE *fplog;
char ofile[MAXPATHLEN];
int i,j,xi,eta,index;
int ix,iy,ir,iz,ntot;
int data_nx,data_ny,data_pl;
char tmpstring[MAXPATHLEN], lognm[MAXPATHLEN];
unsigned nn[2];
float fval,deltaxy,deltaxy1,deltar,r_hat,r,x,y;
float *realNN,*imagNN;
float temp;
float *rowData2;
float *rowData;

ratfactor = 1.0 - 2.0 * AmpRatio;
fprintf(stderr,"ratio factor = %f\n",ratfactor);
imagNN = (float*)head_re->data;
realNN = (float*)head_im->data;

data_nx = (int)(head_im->nx);
data_ny = (int)(head_im->ny);
data_pl = data_nx * data_ny;
Nz = (int)(head_im->nz); 
N = data_nx;
deltaxy = sampling_xy;
fmax = 1/ deltaxy ;
fmax_half = fmax / 2.;

/* because I want to use the same factor for the DIC as in the 64x64 case I have
to compensate for the fact the image is 128x128 that I need for the 
rotation. So it seems that I need to multiply by a factor of 2 */ 
deltaf = 2. * fmax / N;
delta = 2. / N;

fprintf(stderr,"deltaxy= %f (mm), deltaf= %f (1/mm)\n", deltaxy,deltaf);

if((fplog=fopen("log.info","w"))==(FILE*)NULL)
  {
  fprintf(stderr,"ERROR: Can't open '%s' for write/append\n",lognm);
  exit(1);
  }
fprintf(fplog,"                 File to save PSF : %s\n",ofile);
fprintf(fplog,"Shear along x in image space (mm) : %f\n",DeltaX);
fprintf(fplog,"                Phase bias (rads) : %f\n",DeltaPhi);
fprintf(fplog,"                  Amplitude ratio : %f\n",AmpRatio);
fprintf(fplog,"   Pixel size in image space (mm) : %f\n",deltaxy);
fprintf(fplog," Pixel size in freq. space (1/mm) : %f\n",deltaf);
fprintf(fplog,"\n");
fclose(fplog);
DeltaPhi = DeltaPhi / 2.0;
DeltaX = DeltaX / 2.0;

num_pix = (int)(DeltaX / deltaxy);

twoPInum_pix = 2.*PI*num_pix;

fprintf(stderr,"DeltaX/2= %f mm (num_pix=%d) DeltaPhi/2=%f rads\n",DeltaX, num_pix ,DeltaPhi);

printf("twoPInum_pix=%f, DeltaPhi=%f\n",twoPInum_pix,DeltaPhi);
fflush(stdout);

if (ratfactor == 0.0) {
  for(i=0;i<N;i++)
  {
  ir= i * N;
  for  (j=0;j<N;j++)
   {
      
      index=ir + j;

      f = (float) j * deltaf; 

      if (f > fmax_half)   
         f -= fmax;

      f = f/fmax;

      /* multiply by phase factor for DIC optics */
      /* phase factor was changed on 1/31/95; also on 7/30/96 I added the "-"
         do be consistent with the FFT exp sign. */

      sinfactor =  sin(twoPInum_pix * f + DeltaPhi) ;
      
      temp =  sinfactor * outpsf_im[index];
      outpsf_im[index] = - sinfactor * outpsf_re[index];
      outpsf_re[index] = temp;
    }
   }
 } else {
 	for(i=0;i<N;i++)
	  {
            ir= i * N;
  	    for(j=0;j<N;j++)
   	    {
   	    index=ir + j;
   	    f = (float) j * deltaf; 
   	    if (f > fmax_half)	f -= fmax;
	    f = f/ fmax;
   	    /* multiply by phase factor for DIC optics */
   	    angle = twoPInum_pix * f + DeltaPhi;

	    /* note that the rest has been changed as of 3/8/95 */
   	    sinfactor = sin(angle);
   	    cosfactor = cos(angle) * ratfactor;
   	    temp = cosfactor * outpsf_re[index] + sinfactor * outpsf_im[index];
   	    outpsf_im[index] = cosfactor*outpsf_im[index] - sinfactor*outpsf_re[index] ;
   	    outpsf_re[index] = temp;
    	    }
       }
    }

printf("index is %d,N is %d!!\n",index,N);
printf("Leaving!!\n");
fprintf(stderr,"DONE.\n");
fflush(stdout);
}
