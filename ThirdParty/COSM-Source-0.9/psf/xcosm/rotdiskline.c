/***************************************************************************
  COPYRIGHT 2000 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
/**************************************************************

rotdiskline.c

This program takes an X-Z crossection of the 3-D intensity 
point-spread function of the optical microscope and calculates 
a full 3-D (XYZ) psf for the rotating disk partially confocal 
microscope (PCSM) with slit appertures.  It is assumed that the 
condenser/objective lens PSF has circular symmetry. And that the slits
are horiontal (parallel to the x axis)

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

extern void WritePlane(int plane, float *outvol, float *inslice, int w, int h, int inx, int iny);
extern void ZeroOut(float *f, int siz);

extern float deltar; /* real deltar value from calling routine */

void ScaleVol(float *f,int siz, float fsc)
{
int i;
float *fptr;

if(fsc == 0.0) fsc = 1.0;
for(fptr=f,i=0;i<siz;i++,fptr++) {
   *fptr /= fsc;
//printf("ScaleVol; i: %d, val: %d\n", i, *fptr);
}

}

void rotdiskline(float *outpsf, osm_ds *head, int bin, 
		float distance, float size, int tandem)
{
int	Nxy,Nz,Nr;
int	ovrsmp;
int	symmetric;
FILE 	*fpp,*logfp,*infofp;

float   *psfxz;
float 	*psf;		/* One (xy) plane of the illumination PSF              */
float 	*psfbin; 	/* One plane of the condenser PSF and of the final PSF */
float 	*psfcond;	
float   *slicerow; 	/* 1D array used for cubic spine interpolation */
double  *bcf, *ccf, *dcf, *dr_cord; /* arraies for cubic spline inter. */
double  cord, sq_cord;
int endpoints;
float interp_val;

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
int 	Mxy, Halfway;
int 	Nlsf; 		/* Number of samples in the line spread before interp  */
int 	Mlsf;
float 	nrm_deltar, 	/* radial or lateral sampling rate                     */
	deltaxy=1.0, deltaxy1;
float 	weight; 	/* Weigth factor for apodization                       */
float 	r, rd, x, y, xsq; 	/* Distances in the detector plane (mm)                */
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
deltaxy   = head->ylength;
*/

if(distance <= 0.0) distance = 1E20;

/* Initialize a couple fo things that depend on the prompted data */
HalfNxy = Nxy/2;

/* Number of pixels before binning */
Mxy = Nxy*bin;
Halfway = (HalfNxy+1)*bin;

/* Number of points in the lsf before interpolation to a coarser grid */
Nlsf = (int)( ((float)(Nr-1))/sqrt(2.0));

/* Minimun amd maximum values for index in summation of line spreads */
Minn1 = (int)( ((float)(-Mxy))/distance);
Maxn1 = (int)( ((float)(2*Mxy))/distance);

/* radial increment */
nrm_deltar = deltaxy/((float)ovrsmp);

xmax = ((float)(Halfway-1))*deltaxy;
Maxxprime = amax1((double)( ((float)Maxn1)*distance), 
                   amax1(myfabs(xmax-((float)Minn1)*distance), 
                         myfabs(xmax-((float)Maxn1)*distance)));
Maxir = Maxxprime/nrm_deltar + 1;

/* Scale distance between apertures by pre-binning pixel size */
distance = distance*deltaxy;

Zup = Nz;
if (symmetric) Zup = (Nz/2) + 1;

/*   ..Make sure that the log file exists */

if((logfp=fopen(lognm,"a"))==(FILE *)NULL) {
   fprintf(stderr,"WARNING: can't open logfile %s for write.\n",lognm);
   logfp = stderr;
   }
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


fprintf(stderr,"check file %s for progress report.\n",lognm);
     
/*  Calculate smallest power of two grater than or equal to 2*Nlsf */
ix = intlog2(2*Nlsf-1); 
Mlsf=(int)(pow(2.0,(double)(ix+1)));

/* allocate memory for slit array */
if((slit=(float *)calloc(sizeof(float),(Mlsf+2)))==(float *)NULL){
   fprintf(stderr,"out of memory in %s\n",prognm);
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
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(work,Nlsf);

if((radpsf=(float *)calloc(sizeof(float),Nr))==(float *)NULL){
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(radpsf,Nr);

if((lsf=(float *)calloc(sizeof(float),(Mlsf+2)))==(float *)NULL){
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(lsf,Mlsf+2);
Clsf = (fcomplex *)lsf;

if((psfcond=(float *)calloc(sizeof(float),Mxy*Mxy))==(float *)NULL){
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(psfcond,Mxy*Mxy);
 
if((psf=(float *)calloc(sizeof(float),Mxy))==(float *)NULL){
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(psf,Mxy);

if((psfbin=(float *)calloc(sizeof(float),Nxy*Nxy))==(float *)NULL){
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(psfbin,Nxy*Nxy);


/* allocate arrays for the cubic spline interpolation routine (spline.c) */

if((dr_cord=(double*)calloc(sizeof(double),Nr))==(double*)NULL) {
 fprintf(stderr,"ERROR: Out of memory (6).\n");
 exit(1);
 }
if((bcf=(double*)calloc(sizeof(double),Nr))==(double*)NULL) {
 fprintf(stderr,"ERROR: Out of memory (6).\n");
 exit(1);
 }
if((ccf=(double*)calloc(sizeof(double),Nr))==(double*)NULL) {
 fprintf(stderr,"ERROR: Out of memory (6).\n");
 exit(1);
 }
if((dcf=(double*)calloc(sizeof(double),Nr))==(double*)NULL) {
 fprintf(stderr,"ERROR: Out of memory (6).\n");
 exit(1);
 }

slicerow = (float*) psfxz;
/* initialize radial coordinates array for the spline call */
 for(ir=0;ir<Nr;ir++) 
	dr_cord[ir] = ir * deltar; /* use real deltar for the interpolation */

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

   for(ix=2;ix<=Nlsf;ix++)
      *(lsf+Mlsf-ix+1) = *(lsf+ix-1);

   real_fft3d(lsf,Mlsf,1,1,Mlsf+2,1,FORWARD);
   ScaleVol(lsf,Mlsf+2,(float)Mlsf);

   mult3dcm(Clsf,Cslit,Mlsf/2+1);

   real_fft3d(lsf,Mlsf,1,1,Mlsf+2,1,REVERSE);
   for(ix=0;ix<Mlsf+2;ix++) *(lsf+ix) *= (float)Mlsf;


/* compute the cubic spline interpolation coefficients: bcf, ccf, dcf */
spline_(&Nr,dr_cord, slicerow,bcf,ccf,dcf);
deltaxy1 = deltar*(float)ovrsmp;

   /* Interpolate between tabulated values */
   for(ix=1;ix<=Halfway;ix++){
      x = ((float)(ix-1))*deltaxy1;
      xsq = x*x;
      for(iy=1;iy<=Halfway;iy++){
         y = ((float)(iy-1))*deltaxy1;
	 r = (float)sqrt((double)(xsq + y*y));
	 rd = r/deltar;
	 ir = (int) rd;
	 cord = r - dr_cord[ir];
	 sq_cord = cord * cord;

         /* interpolated value based on cubic spline interpolation equation */
         interp_val = slicerow[ir] + bcf[ir]*cord + ccf[ir]*sq_cord + dcf[ir]*cord*sq_cord;
         *(psfcond+ix-1+(iy-1)*Mxy) =  interp_val;
      }

      /* Initialize to use as accumulator */
      illum = 0.0;
      /* accumulate contributions of different slits on the disk */
      for(n1=Minn1;n1<=Maxn1;n1++){
         xprime = x - distance * ((float)n1);
         rprime = myfabs(xprime)/nrm_deltar;
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

      }

      slicerow += Nr;
   /*      
      Circular replication into the other three quadrants
      NOTE: samples at ix or iy) = HalfNxy+1 do not have replicas
   */
   /* Replicate I into II */
   /* Replicate I and II into IV and III, respectivelly */

   for(ix=2;ix<=bin*HalfNxy;ix++){
      *(psf+Mxy-ix+1) = *(psf+ix-1);
      for(iy=1;iy<=Halfway;iy++)
         *(psfcond+Mxy-ix+1+(iy-1)*Mxy) = *(psfcond+ix-1+(iy-1)*Mxy);
   }
   for(ix=1;ix<=Mxy;ix++){
      for(iy=2;iy<=bin*HalfNxy;iy++)
         *(psfcond+ix-1+(Mxy-iy+1)*Mxy) = *(psfcond+ix-1+(iy-1)*Mxy);
      }

   /* 
      Multiply illumination and condenser PSFs (horizontal slit assumed)
   */
   for(iy=1;iy<=Mxy;iy++) {
      ftmp = *(psf+iy-1);
      for(ix=1;ix<=Mxy;ix++)
         *(psfcond+ix-1+(iy-1)*Mxy) *= ftmp;
      }

   /* _Bin here */
   /* 
      If binning by a factor of 1 the folowing loops 
      copy psfcond to psfbin
   */
   for(ix=1;ix<=Nxy;ix++)
      for(iy=1;iy<=Nxy;iy++)
         *(psfbin+ix-1+(iy-1)*Nxy) = 0.0;
   for(ix=1;ix<=Mxy;ix++){
      jx = ((ix-1)/bin) + 1;
      for(iy=1;iy<=Mxy;iy++){      /* BUG1 - was <=Nxy */
         jy = ((iy-1)/bin) + 1;
         *(psfbin+jx-1+(jy-1)*Nxy) += *(psfcond+ix-1+(iy-1)*Mxy);
         }
      }

   for(iy=1;iy<=Nxy;iy++)
    for(ix=1;ix<=Nxy;ix++)
       if( *(psfbin + ix-1 + (iy-1)*Nxy) < 0.0) 
           *(psfbin + ix-1 + (iy-1)*Nxy) = 0.0;

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
