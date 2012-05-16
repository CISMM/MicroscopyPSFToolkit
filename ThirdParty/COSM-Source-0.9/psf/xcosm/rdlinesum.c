/***************************************************************************
  COPYRIGHT 2000 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
/**************************************************************

rdlinesum.c

This program takes an X-Z crossection of the 3-D intensity 
point-spread function of the optical microscope and calculates 
a full 3-D (XYZ) psf for the rotating disk partially confocal 
microscope (PCSM) with slit appertures.  It is assumed that the 
condenser/objective lens PSF has circular symmetry. And that the slits
are horiontal (parallel to the x axis)
It is used instead of rotdiscline.c when the desirable pixel size in the 
final PSF is larger than the one computed based on the Nyquist rate, i.e.
when the final PSF is undersampled. The routine computes and
oversampled PSF and at the end it sums neighboring pixels (at least 5) to
get the final PSF.


CAVEATS:
The XZ crossection needs a finer sampling rate than the resulting
3-D PSF because the latter is calculated by linearly interpolating
between the samples of the XZ crossection.
The XZ crossection must have wider support radially than the 
PSF of the rotating disk micorscope.  The 3-D PSF of the rotating-disk
PCSM, hp(x,y,z), is the superposition of shifted versions of the 
circularly symmetric PSF h(r,z).  If hc(r,z) has support over a 
region about the size of the region of support of the PSCM's PSF,
there will be artifacts at the edges of hp(x,y,z). A good rule of 
thumb is to give hc(r,z) a region of support 2 to 3 times the 
support desired for hp(x,y,z)
**********************************************************/

#include "misc.h"
#include "washu.h"

extern char prognm[];
extern char psfnm[];
extern char lognm[];
extern float deltar, deltaxy, deltaxy_nyq;
extern void Sum4N4NToNN(float *sorc, float *desti, int Ndim, int osmp);

extern void WritePlane(int plane, float *outvol, float *inslice, int w, int h, int inx, int iny);
extern void ZeroOut(float *f, int siz);
extern void Save2D(float *datt, int NN, char *ffname);
extern void ScaleVol(float *f,int siz, float fsc);

void rdlinesum(float *outpsf, osm_ds *head, int bin, 
		float distance, float size, int tandem)
{
int	Nxy,Nz,Nr;
int	ovrsmp;
int	symmetric;
FILE 	*fpp,*logfp,*infofp;

float   *psfxz;
float 	*psf;		/* One (xy) plane of the illumination PSF              */
float 	*psfbin1, *psfbin; /* One plane of the final and of the condenser PSF  */
float 	*psfcond;		
float 	illum; 		/* One sample of the illumination patern               */
float 	*lsf; 		/* Line spread Function                                */
fcomplex *Clsf;
float 	*slit; 		/* Slit function                                       */
fcomplex *Cslit;
float 	peak; 		/* to normalize                                        */
float 	*work; 		/*  Working storage                                    */
float 	*radpsf; 	/*  PSF intensity image crossection PSF (x, 0, z)      */
int 	iz, ir, ix, 	/* indices for depth, radial distance,                 */
	iy, jx, jy; 	/* and the two lateral coordinates                     */
int 	n1, n2; 	/* Indices for periodic sampling used in the PSF       */
int 	Maxn1, Maxn2, 
	Minn1, Minn2;
int 	HalfNxy, Zup;
int	Lxy;		/* Number of samples in (x,y) before binning  */
int	middle;		/* middle of oversampled array */

/* Number of samples in (x,y) before binning and after 
   oversampling for the summation*/
int 	Mxy, Halfway;   

int 	Nlsf; 		/* Number of samples in the line spread before interp  */
int 	Mlsf;
float 	nrm_deltar, 	/* normalized radial or lateral sampling rate          */
	nrm_deltaxy=1.0;
/*   new lateral sampling rate based on oversampling */
float           deltaxy1;
float 	weight; 	/* Weigth factor for apodization                       */
float 	rd, x, y, xsq; 	/* Distances in the detector plane (mm)                */
/*   oversampling rate for the xy plane */
int             osamp;

float 	xmax;
float 	xprime, yprime, rprime;
float 	Maxxprime;
float   ftmp;
int 	Maxir;
char 	temp_string[255];
char	answer[255];

/* --------- BEGIN EXECUTABLE CODE -------------------------- */

Nr        = head->nx;
Nz        = head->ny;
psfxz     = (float*)head->data;
ovrsmp    = head->xstart;
Nxy       = head->ystart;
symmetric = head->zstart;
/*
nrm_deltar    = head->xlength;
nrm_deltaxy   = head->ylength;
*/

if(distance <= 0.0) distance = 1E20;

/* Initialize a couple fo things that depend on the prompted data */
HalfNxy = Nxy/2;

/* Number of pixels before binning */
Lxy = Nxy*bin;

/* oversampling rate    */
/* make it beat nyquist */
osamp = (int)(deltaxy/deltaxy_nyq) + 1;

if(osamp < 5) osamp = 5;

/* make odd */
if(!(osamp & 1)) osamp++;

/* New deltaxy for oversampling */
deltaxy1 = deltaxy/((float)osamp);

printf("rdlinesum; DXY: %f DXY_NYQ: %f SAMP = %d NEW_DXY: %f\n",
		deltaxy,deltaxy_nyq,osamp,deltaxy1);

Mxy = Lxy * osamp; /* number of pixels after oversampling */
Halfway = (HalfNxy+1)*bin*osamp;
printf("rdlinesum; Mxy = %d \n", Mxy);


/* Number of points in the lsf before interpolation to a coarser grid */
Nlsf = (int)( ((float)(Nr-1))/sqrt(2.0));

/* Minimun amd maximum values for index in summation of line spreads */
Minn1 = (int)( ((float)(-Mxy))/distance);
Maxn1 = (int)( ((float)(2*Mxy))/distance);

/* radial increment */
nrm_deltar = nrm_deltaxy/((float)ovrsmp);

xmax = ((float)(Halfway-1))*nrm_deltaxy;
Maxxprime = amax1((double)( ((float)Maxn1)*distance), 
                   amax1(myfabs(xmax-((float)Minn1)*distance), 
                         myfabs(xmax-((float)Maxn1)*distance)));
Maxir = Maxxprime/nrm_deltar + 1;

/* Scale distance between apertures by pre-binning pixel size */
distance = distance*nrm_deltaxy;

Zup = Nz;
if (symmetric) Zup = (Nz/2) + 1;

/*   ..Make sure that the log file exists */

if((logfp=fopen(lognm,"a"))==(FILE *)NULL) {
   fprintf(stderr,"rdlinesum: WARNING: can't open logfile %s for write.\n",lognm);
   logfp = stderr;
   }
fprintf(logfp," rdlinesum;\n");
fprintf(logfp,"         Number of slices or planes: %d\n", Nz);
fprintf(logfp," Number of xy-samples after binning: %d\n", Nxy);
fprintf(logfp," Number of radial samples available: %d\n", Nr);
fprintf(logfp,"   Oversampling (before binning) by: %d\n", ovrsmp);
fprintf(logfp,"                     binning factor: %d\n", bin);
fprintf(logfp,"    distance between slit apertures: %E\n", distance);
fprintf(logfp,"width of slit apertures (ovrsmpled): %f\n", size);
fprintf(logfp,"even symmetry in z is %s assumed.\n",
  symmetric ? " ":"NOT");

/*   ..write info file for PSF */
   sprintf(temp_string,"%s.info",psfnm);
   if((infofp=fopen(temp_string,"a"))==(FILE *)NULL) 
      fprintf(stderr,"can't open file %s for write.\n",temp_string);
   else {
      sprintf(temp_string,"new psf file: %s",psfnm);
      keeplog (infofp, temp_string);
      fprintf(infofp,"        Number of slices or planes: %d\n", Nz);
      fprintf(infofp,"Number of xy-samples after binning: %d\n", Nxy);
      fprintf(infofp,"Number of radial samples available: %d\n", Nr);
      fprintf(infofp,"  Oversampling (before binning) by: %d\n", ovrsmp);
      fprintf(infofp,"                    binning factor: %d\n", bin);
      fprintf(infofp,"   distance between slit apertures: %E\n", distance);
      fprintf(infofp,"   aperture width(prebin,oversamp): %f\n", size);
      fprintf(infofp,"even symmetry in z is %sassumed.\n",
   	   symmetric ? " ":"NOT ");
      fclose(infofp);
      }


fprintf(stderr,"rdlinesum; heck file %s for progress report.\n",lognm);
     
/*  Calculate smallest power of two grater than or equal to 2*Nlsf */
ix = intlog2(2*Nlsf-1); 
Mlsf=(int)(pow(2.0,(double)(ix+1)));

/* allocate memory for slit array */
if((slit=(float *)calloc(sizeof(float),(Mlsf+2)))==(float *)NULL){
   fprintf(stderr,"rdlinesum; out of memory in %s\n",prognm);
   exit(1);
}
ZeroOut(slit,Mlsf+2);

/* Sample the slit funtion */
for(ix=0;ix<Mlsf;ix++) *(slit+ix) = 0.0;
*slit = 1.0;
for(ix=1;ix<=( (int)(size/2.0) );ix++){
   *(slit+ix) = 1.0;
   *(slit+Mlsf-ix) = 1.0;
    }

/* Correction for even number of pixels */
if ( ((int)size)%2 == 0) *(slit+Mlsf- (int)(size/2.0)) = 0.0;

/* Fourier transform the slit function */
real_fft3d(slit, Mlsf,1,1,Mlsf+2,1,FORWARD);
ScaleVol(slit,Mlsf+2,(float)Mlsf);

Cslit = (fcomplex *)slit;
peak = 0;

/* 
   If tandem scanning, the equivalent slit function 
   is the convolution of the slit with itself
*/
if (tandem) mult3dcm (Cslit, Cslit, Mlsf/2+1);
/* since we are multiplying the complex array Cslit
with itself, there are only Mlsf/2+1 elements,
instead of Mlsf+2 elements for the real array slit. */

/* 
   allocate memory for some arrays
   before entering iz loop 
*/
if((work=(float *)calloc(sizeof(float),Nlsf))==(float *)NULL){
   fprintf(stderr,"rdlinesum; out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(work,Nlsf);

if((radpsf=(float *)calloc(sizeof(float),Nr))==(float *)NULL){
   fprintf(stderr,"rdlinesum; out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(radpsf,Nr);

if((lsf=(float *)calloc(sizeof(float),(Mlsf+2)))==(float *)NULL){
   fprintf(stderr,"rdlinesum; out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(lsf,Mlsf+2);
Clsf = (fcomplex *)lsf;

if((psfcond=(float *)calloc(sizeof(float),Mxy*Mxy))==(float *)NULL){
   fprintf(stderr,"rdlinesum; out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(psfcond,Mxy*Mxy);
 
if((psf=(float *)calloc(sizeof(float),Mxy))==(float *)NULL){
   fprintf(stderr,"rdlinesum; out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(psf,Mxy);

if((psfbin=(float *)calloc(sizeof(float),Nxy*Nxy))==(float *)NULL){
   fprintf(stderr, "rdlinesum; out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(psfbin,Nxy*Nxy);

for(iz=1;iz<=Zup;iz++){
   sprintf(temp_string,"iz = %d",iz);
   keeplog(logfp, temp_string);

   radpsf = (psfxz + (iz-1)*Nr);

   /* 
      The following loops are based on the 
      PSF's even symmetry in x and y
   */
   for(ix=1;ix<=Nlsf;ix++){
      /* 
         x, y, and r have the same resolution in this loop. 
         No scaling needed
      */
      x = ix-1;
      xsq = x*x;
      for(iy=1;iy<=Nlsf;iy++){
         /* 
           ..LSF is sampled at the same resolution as the radial cut
         */
         y = iy-1;
	 rd = (float)sqrt((double)(xsq + y*y));
         ir = rd+1;

         /*   ..Interpolate */
	 *(work+iy-1) = *(radpsf+ir-1) * (((float)ir) - rd) 
                        + *(radpsf+ir) * (rd + 1.0 - ((float)ir));
         }
      /*  __Integrate PSF to get Line Spread Function (lsf) */
      /* Initialize trapezoidal rule */
      *(lsf+ix-1) = - 0.5 * ( (*work) + *(work+Nlsf-1)); 

      /* Accumulate  */
      for(iy=1;iy<=Nlsf;iy++)
	       *(lsf+ix-1) += *(work+iy-1);

      /* Scale by the integration step (equal to the radial sampling)*/
      *(lsf+ix-1) *= nrm_deltar;
   }

   for(ix=2;ix<=Nlsf;ix++) {
      *(lsf+Mlsf-ix+1) = *(lsf+ix-1);
   }

   real_fft3d(lsf,Mlsf,1,1,Mlsf+2,1,FORWARD);
   printf("real_fft3d\n");
   ScaleVol(lsf,Mlsf+2,(float)Mlsf);

   mult3dcm(Clsf,Cslit,Mlsf/2+1);
   printf("multi3dcm\n");
   real_fft3d(lsf,Mlsf,1,1,Mlsf+2,1,REVERSE);
   for(ix=0;ix<Mlsf+2;ix++) {
	*(lsf+ix) *= (float)Mlsf;
   }

   printf("rdlinesum: deltar: %g, deltaxy1: %g \n", deltar, deltaxy1);

   /* Interpolate between tabulated values */
   for(ix=1;ix<=Halfway;ix++){
      x = ((float)(ix-1))*deltaxy1;
      xsq = x*x;
      for(iy=1;iy<=Halfway;iy++){
         y = ((float)(iy-1))*deltaxy1;
	 rd = (float)sqrt((double)(xsq + y*y))/deltar;
	 ir = rd+1;
         *(psfcond+ix-1+(iy-1)*Mxy) =  *(radpsf+ir-1) * (((float)ir) - rd)
                            + *(radpsf+ir) * (rd + 1.0 - ((float)ir));
         } /* end loop for iy */

      /* Initialize to use as accumulator */
      illum = 0.0;
      /* accumulate contributions of different slits on the disk */
      for(n1=Minn1;n1<=Maxn1;n1++){
         xprime = x - distance * ((float)n1);
	 /* convert resolution from radial sampling to final xy sampling */
         rprime = myfabs(xprime) * (deltaxy1 / deltar);
         ir = rprime + 1;

         /* Interpolate and accumulate */
         if (ir<Nlsf) {
            /* Interpolate between adjacent samples */
            illum =illum + *(lsf+ir-1) * (((float)ir)-rprime)
                         + *(lsf+ir)   * (rprime+1.0-((float)ir));
            }
         else {
            /* The following if () should be dropped at some point */
	    if (ir>Maxir) {
	       ir = Maxir;
	       printf("ir was set = %d\n", Maxir);
	       printf("Check code again\n");
               }
            /* 
               Apodize from last available sample to zero
               to avoid sharp edges in the PSF
               Apodize with a linear interpolation to zero at the 
               largest value of ir expected.
            */
            weight = (float)(Maxir-ir)/(float)(Maxir-Nlsf);
            illum += ( *(lsf+Nlsf-1) * weight );
            }

         }

      if (illum<0.0)
         printf("illum<0 at (ix,iy,iz,n1)=%d %d %d %d\n",ix,iy,iz,n1);

      *(psf+ix-1) = illum ;

      } /* end loop for ix */

   /*      
      Circular replication into the other three quadrants
      NOTE: samples at ix or iy) = HalfNxy+1 do not have replicas
   */
   /* Replicate I into II */
   /* Replicate I and II into IV and III, respectivelly */
   middle = bin*HalfNxy*osamp;
   for(ix=2;ix<=middle;ix++){
      *(psf+Mxy-ix+1) = *(psf+ix-1);
      for(iy=1;iy<=middle+1;iy++)
         *(psfcond+Mxy-ix+1+(iy-1)*Mxy) = *(psfcond+ix-1+(iy-1)*Mxy);
   }
   for(ix=1;ix<=Mxy;ix++){
      for(iy=2;iy<=middle;iy++)
         *(psfcond+ix-1+(Mxy-iy+1)*Mxy) = *(psfcond+ix-1+(iy-1)*Mxy);
      }

   /* 
      Multiply illumination and condenser PSFs (horizontal slit assumed)
   */
   for(iy=1;iy<=Mxy;iy++) {
      ftmp = *(psf+iy-1);
      for(ix=1;ix<=Mxy;ix++) {
         *(psfcond+ix-1+(iy-1)*Mxy) *= ftmp;
      }
   }

   /* allocate memory for psfbin1. Notice size of array */
	if((psfbin1=(float *)calloc(sizeof(float),Lxy*Lxy))==(float *)NULL)
	{
	   fprintf(stderr,"rdlinesum; out of memory in %s\n",prognm);
	   exit(1);
	 }

    /* sum neighboring pixels to downsample the BIG PSF */
    /* note dimension goes from Mxy to Lxy = Nxy*bin */

    Sum4N4NToNN(psfcond,psfbin1,Lxy,osamp); /* routine is in rotsum.c */

   /* _Bin here */
   /* 
      If binning by a factor of 1 the folowing loops 
      copy psfbin1 to psfbin
   */
   if (bin == 1) /* Lxy = Nxy */
	cp3dr (psfbin, psfbin1, Nxy*Nxy);
   else {
	 /* Add pixels together */
         /* Clean psfbin array to use as accumulator */
   	for(ix=1;ix<=Nxy;ix++)
	      for(iy=1;iy<=Nxy;iy++)
         	*(psfbin+ix-1+(iy-1)*Nxy) = 0.0;
	 /* binning */
   	for(ix=1;ix<=Lxy;ix++){
	      jx = ((ix-1)/bin) + 1;
	      for(iy=1;iy<=Lxy;iy++){      
	         jy = ((iy-1)/bin) + 1;
	         *(psfbin+jx-1+(jy-1)*Nxy) += *(psfbin1+ix-1+(iy-1)*Lxy);
              }
         }
   }
   for(iy=1;iy<=Nxy;iy++)
    for(ix=1;ix<=Nxy;ix++) {
       if( *(psfbin + ix-1 + (iy-1)*Nxy) < 0.0) 
           *(psfbin + ix-1 + (iy-1)*Nxy) = 0.0; 
    }

   WritePlane(iz-1,outpsf,psfbin,Nxy,Nxy,Nxy,Nxy);
   if(symmetric && (iz > 1))
       WritePlane(Nz-iz+1,outpsf,psfbin,Nxy,Nxy,Nxy,Nxy);


   }  /* end of "for(iz=1;iz<=Zup;iz++)" */


/* calculate and save line spread function if desired */

keeplog (logfp, "DONE");
if(logfp != stderr) fclose(logfp);

/*
free(radpsf);
free(lsf);
free(psfbin);
free(psfcond);
free(psf);
free(work);
free(slit);
*/
}
