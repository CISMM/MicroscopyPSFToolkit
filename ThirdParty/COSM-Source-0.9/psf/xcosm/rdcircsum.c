/***************************************************************************
  COPYRIGHT 2000 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
/***********************************************************

rdcircsum.c

This program takes an X-Z crossection of the 3-D intensity 
point-spread function of the non-confocal optical microscope and
calculates a 3-D (XYZ) psf for the rotating disk partially confocal 
microscope (PCSM) with finite sized pinhole apertures.
It is used instead of rotdisccirc.c when the desirable pixel size in the
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
***********************************************************/

#include "misc.h"
#include "washu.h"

extern float deltaxy,deltaxy_nyq;
extern float deltar;
extern char prognm[];
extern char psfnm[];
extern char lognm[];
extern void WritePlane(int plane, float *outvol, float *inslice, int w, int h, int inx, int iny);
extern void ZeroOut(float *f,int siz);
extern void Sum4N4NToNN(float *sorc, float *desti, int Ndim, int osmp); 

void rdcircsum(float *outpsf, osm_ds *head, int bin,
		float distance, float size, int tandem)
{
int 		Nxy,Nz,Nr;
int		ovrsmp;
int		symmetric;
FILE		*logfp,*infofp;

float		*pupil; 	/*   Pupil function and its Fourier transform */
fcomplex	*Cpupil;
float		*psfobj; 	/*   One (xy) plane of the objective's PSF    */
float		*psfbin1; 	/*   One plane (xy) of the intermediate PSF   */
float		*psfbin; 	/*   One plane (xy) of the final PSF          */
float		*psfcond; 	/*   One plane of the condenser PSF           */
fcomplex	*otfcond;
float		peak;		/*   Max. value of arrays (normalization)     */
float		illum; 		/*   One sample of the illumination patern    */
float		tmp;
float		alpha;
float		*radpsf; 	/*   PSF intensity image XY crossection       */
float		*psfxz;
float		Vxx, Vyx, 	/*   Elements of the sampling matrix          */
		Vxy, Vyy; 
int		iz, ir, ix, 	/*   indices for depth, radial distance,      */
		iy, jx, jy; 	/*    and the two lateral coordinates 	      */
int		n1, n2; 	/*   indices for periodic sampling in the PSF */
int		Maxn1, Maxn2, 
		Minn1, Minn2;
int		Maxir;

/* __Stuff to calculte min and max index values */
float		Sxmin, Sxmax, xmax, Maxxprime;
float		Symin, Symax, ymax, Maxyprime;
float		Maxrprime;
int		sMinn1, lMaxn1;
float		mytemp;

/*   Number of samples in (x, y), in z and radially (= Nxy*sqrt(2))*/
int		HalfNxy, Zup;
/*   Number of samples radially on the convolution of psfcond and pupil */
int		NetNr;
/*   Number of samples in (x,y) before binning  */
int		Lxy;
/*   Number of samples in (x,y) before binning and after 
     oversampling for the summation*/
int		Mxy, Halfway;
/*   Pupil function region of support */
int		Pxy;
/*   normalized radial or lateral sampling rate */
float		nrm_deltar, nrm_deltaxy=1.0;
/*   new lateral sampling rate based on oversampling */
float		deltaxy1;
/*   oversampling rate for the xy plane */
int		osamp;
/*   Floating point version of the above */
float		ratio;
/*   Weigth factors for apodization */
float		weight, appod;
/*   Distances in the detector plane (mm) */
float		rd, x, y, ysq;
/*   Shifted coordinates  */
float		xprime, yprime, rprime;
/*   Shifting of coodrinates */
float		dx, dy;
/*   Shifted radial coordinate normalized by deltar */
float		normRprime;
    
/*   ..To prompt for things */
char	temp_string[255];
char	answer[255];

/* --------------------- BEGIN EXECUTABLE CODE -------------------- */


/* __Initialize whatever paremeters here */
/*   ..Periodicity matrix */
Vxx = 1.0;
Vyx = 0.0;
Vxy = -0.5;
Vyy = (float)sqrt(3.0)*0.5;

/*   ..Prompt for all data */

Nr	  = head->nx;
Nz	  = head->ny;
psfxz     = (float*)head->data;
ovrsmp    = head->xstart;
Nxy       = head->ystart;
symmetric = head->zstart;
/* leave at 1.0
nrm_deltar    = head->xlength;
nrm_deltaxy   = head->ylength;
*/


HalfNxy = Nxy/2;
Lxy = Nxy*bin; 			/* Number of pixels before binning */
Halfway = (HalfNxy+1)*bin; 
NetNr = Nr; 			/* Number of pixels in convolution of psfcond and pupil */

nrm_deltar = nrm_deltaxy/((float)ovrsmp);

/* oversampling rate    */
/* make it beat nyquist */
osamp = (int)(deltaxy/deltaxy_nyq) + 1;

if(osamp < 5) osamp = 5;

/* make odd */
if(!(osamp & 1)) osamp++;

/* New deltaxy for oversampling */
deltaxy1 = deltaxy/((float)osamp);
ratio = (float) deltaxy1 / deltar;

printf("DXY: %.4f DXY_NYQ: %.4f SAMP = %d NEW_DXY: %.4f\n",
		deltaxy,deltaxy_nyq,osamp,deltaxy1);

Mxy = Lxy * osamp; /* number of pixels after oversampling */
Halfway = (HalfNxy+1)*bin*osamp; 
printf("ratio = %f Mxy = %d \n", ratio,Mxy);

/* check for single aperture */
if(distance <= 0.0) distance = 1E20;

distance = distance*nrm_deltaxy; 	/* Scale by pre-binning pixel size */

Zup = Nz; 			/* Limit for z index */
if (symmetric) Zup = (Nz/2) + 1;

/*   ..Maximum values for row sampling index 
	(NOTE: Assumes Vyy != 0 <--) */
Minn2 = ((float)-Mxy)/(Vyy*distance);
Maxn2 = ((float)(2*Mxy))/(Vyy*distance);

/* __Calculate a conservative extimate of the maximum value if the 'ir'
index that is expected to be used.  The estimate is calculated
with the following assumptions:
Vxx > 0, Vyy > 0 (strictly greater)
Vxy < 0 	 (Striclty less than)
Vyx = 0		 (exactly)
..Most negative value of n1 expected */

sMinn1 = ((float)-Mxy)/distance - (float)Minn2*Vxy/Vxx;
/*   ..Most positive value of n1 expected  */
lMaxn1 =  (((float)(2*Mxy))/distance - (float)Maxn2*Vxy)/Vxx;
/*   ..most negative value of term subtracted from xprime */
Sxmin = Vxx*(float)sMinn1 + Vyx*(float)Maxn2;
/*   ..Most positive value of term subtracted from xprime */
Sxmax = Vxx*(float)lMaxn1 + Vyx*(float)Minn2;
/*   ..Largest value of x */
xmax = (float)(Halfway-1) * nrm_deltaxy;
/*   ..Largest absolute value of xprime  */
Maxxprime = amax1(Sxmax*distance, xmax - Sxmin*distance);
/*   ..Most negative value of the term subtracted 
	from yprime (with Vyx=0,Vyy>0) */
Symin = Vyy*(float)Minn2;
/*   ..Most positive value of the term subtracted 
	from yprime (with Vyx=0,Vyy>0) */
/* In the original fortran code, Maxn1 is used in the next line
without being initialized first. I suspect it's a bug, and should
really be Maxn2. */
Symax = Vyy*(float)Maxn2;
/*   ..largest value of y  */
ymax = (float)(Halfway-1) * nrm_deltaxy;
/*   ..Largest absolute value of yprime */
Maxyprime=amax1(Vyy*(float)Maxn2*distance,ymax-Vyy*(float)Minn2*distance);
/*   ..Largest value of rprime we can conservativelly expect */
Maxrprime = (float) sqrt (Maxxprime*Maxxprime + Maxyprime*Maxyprime);
/*   ..Largetst index ir expected */
Maxir = (int) ((Maxrprime + 1.0)/nrm_deltar);
/*   ..Factor to be used in linear interpolation used as apodization */
appod = 1.0;
if (Maxir > Nr) appod = 1.0/(float)(Maxir - Nr) ;

/*   ..Calculate size of pupil function array at the higher resolution*/
Pxy = Nxy*ovrsmp*bin;
/*   ..Calculate the least power of 2 that will fit */
Pxy = intlog2(Pxy-1);
Pxy = (int) pow(2.0 , Pxy+1);

/*   ..Make sure that the log file exists */
if((logfp=fopen(lognm,"a"))==(FILE *)NULL) {
   fprintf(stderr,"WARNING: (rotdiskcirc) can't open logfile `%s' for write.\n",lognm);
   logfp = stderr;
   }
fprintf(logfp,"         Number of slices or planes: %d\n", Nz);
fprintf(logfp," Number of xy-samples after binning: %d\n", Nxy);
fprintf(logfp," Number of radial samples available: %d\n", Nr);
fprintf(logfp,"   Oversampling (before binning) by: %d\n", ovrsmp);
fprintf(logfp,"                     binning factor: %d\n", bin);
fprintf(logfp,"         distance between apertures: %E\n", distance);
fprintf(logfp,"   aperture size(prebin,oversmpled): %f\n", size);
fprintf(logfp,"   pupil support(prebin,oversmpled): %d\n", Pxy);
fprintf(logfp,"                      sample matrix:\n");
fprintf(logfp,"                                Vxx: %f\n",Vxx);
fprintf(logfp,"                                Vxy: %f\n",Vxy);
fprintf(logfp,"                                Vyx: %f\n",Vyx);
fprintf(logfp,"                                Vyy: %f\n",Vyy);
fprintf(logfp,"even symmetry in z is %s assumed.\n", symmetric ? "":"NOT");


   sprintf(temp_string,"%s.info",psfnm);
   if((infofp=fopen(temp_string,"a"))==(FILE *)NULL) 
      fprintf(stderr,"WARNING: can't open info file %s for write.\n",temp_string);
   else {
        fprintf(infofp,"                    binning factor: %d\n", bin);
        fprintf(infofp,"        Distance between apertures: %E\n", distance);
        fprintf(infofp,"aperture diameter(prebin,oversamp): %f\n", size);
        fprintf(infofp,"                     sample matrix:\n");
        fprintf(infofp,"                               Vxx: %f\n",Vxx);
        fprintf(infofp,"                               Vxy: %f\n",Vxy);
        fprintf(infofp,"                               Vyx: %f\n",Vyx);
        fprintf(infofp,"                               Vyy: %f\n",Vyy);
        fclose(infofp);
         }


fprintf(stderr,"check file %s for progress report.\n",lognm);

if (size >= 1.0) {
   keeplog (logfp, "sampling pupil function ");

   /* first allocate memory for array pupil */
   if((pupil=(float *)calloc(sizeof(float),((Pxy+2)*Pxy)))==
      (float *)NULL) {
      fprintf(stderr,"out of memory in %s\n",prognm);
      exit(1);
      }
   ZeroOut(pupil,(Pxy+2)*Pxy);

   /* Get pupil function */
   getpupil(pupil, size, Pxy, Pxy);
   add2columns(pupil,Pxy,Pxy);

   /* Fourier transform pupil function */
   real_fft3d(pupil,Pxy,Pxy,1,Pxy+2,Pxy,FORWARD);
   Cpupil=(fcomplex *)pupil;

   /* If tandem scanning, the equivalent pinhole function is the
      convolution of the pinhole with itself */
   if (tandem) mult3dcm (Cpupil, Cpupil, (Pxy+2)*Pxy/2);

   /* Normalize Cpupil array to avoid handling very small numbers */
   peak = 0.0;
   norm3dcm (Cpupil, &peak, Pxy*(Pxy+2)/2);
   }
     
/* allocate memory for psfobj. Notice size of array   */
if((psfobj=(float *)calloc(sizeof(float),Mxy*Mxy))==(float *)NULL) {
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(psfobj,Mxy*Mxy);

printf("Mxy = %d and Pxy=%d\n",Mxy,Pxy);
fflush(stdout);
/* allocate memory for psfcond. Notice size of array   */
if (Mxy > Pxy){  /* when ovrsmp < 5 then this is true */
  printf("Mxy >Pxy\n");
  if((psfcond=(float *)calloc(sizeof(float),(Mxy+2)*Mxy))==(float *)NULL) {
     fprintf(stderr,"out of memory in %s\n",prognm);
     exit(1);
     }
}
else
  if((psfcond=(float *)calloc(sizeof(float),(Pxy+2)*Pxy))==(float *)NULL) {
     fprintf(stderr,"out of memory in %s\n",prognm);
     exit(1);
     }

ZeroOut(psfcond,(Pxy+2)*Pxy);
otfcond=(fcomplex *)psfcond;

/* allocate memory for psfbin. Notice size of array   */
if((psfbin=(float *)calloc(sizeof(float),Nxy*Nxy))==(float *)NULL) {
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
   }
ZeroOut(psfbin,Nxy*Nxy);

/* allocate memory for psfbi1. Notice size of array */
if((psfbin1=(float *)calloc(sizeof(float),Lxy*Lxy))==(float *)NULL) {
   fprintf(stderr,"out of memory in %s\n",prognm);
   exit(1);
}
ZeroOut(psfbin1,Lxy*Lxy);


for(iz=1;iz<=Zup;iz++){
   if ( (iz-1)%1 == 0) {
      sprintf(temp_string,"iz = %d",iz);
      keeplog(logfp, temp_string);
      }

   /*      ..Read on row of the rz or xz cross-section */
   /* first have to go to the right place in file */

   radpsf = (psfxz + (iz-1)*Nr);

   /* Calculate the objective PSF (iterpolate rz section into xyz sect)*/
   r2xy(psfobj, radpsf, Mxy, Mxy, Nr, ratio);

   if (size >= 1.0) {
      /* Convolve condenser PSF and aperture */
      /* interpolate w/o subsampling */
      r2xy(psfcond, radpsf, Pxy, Pxy, Nr, 1.0);
      add2columns(psfcond,Pxy,Pxy);

      /*  Fourier transform, multiply, inverse Fourier transform */
      real_fft3d(psfcond,Pxy,Pxy,1,Pxy+2,Pxy,FORWARD);
      mult3dcm(otfcond, Cpupil, (Pxy+2)*Pxy/2);
      real_fft3d(psfcond,Pxy,Pxy,1,Pxy+2,Pxy,REVERSE);

      /* Copy first line of convolution to radpsf to use for illumination */
      NetNr = Pxy/2+1;
      cp3dr(radpsf, psfcond, NetNr);
      
   }


   /* The following loops are based on the PSF's even symmetry 
       in x and y */
   /* Calculate the illumination distribution */

   for(iy=1;iy<=Halfway;iy++){
      y = (float)(iy-1)*nrm_deltaxy;
      for(ix=1;ix<=Halfway;ix++){
         x = (float)(ix-1)*nrm_deltaxy;
         /*Initialize to use as accumulator */
	 illum = 0.0;
         /*  ..Index n2 is for columns  */
         for(n2=Minn2;n2<=Maxn2;n2++){
            dx = Vxy*n2;
            dy = Vyy*n2;
            
            /* ..These two lines assume Vxx is not zero */
            Minn1=(int)(((float)-Mxy/distance-dx)/Vxx);
            Maxn1=(int)((2.0*(float)Mxy/distance-dx)/Vxx);

            /*  ..Index n1 involves rows and columns */
            for(n1=Minn1;n1<=Maxn1;n1++){
               xprime = x - distance*(Vxx*(float)n1+dx);
               yprime = y - distance*(Vyx*(float)n1+dy);

               /* Rectangular sampling to check what's going on here */
               /*
                          xprime = x - distance*float(n2)
                          yprime = y - distance*float(n1)
               */
               rprime = (float)sqrt((double)(xprime*xprime) + (double)(yprime*yprime));
               /*  ..Normalize to the subsampled grid */
               normRprime = rprime*ratio;
               ir = (int)normRprime + 1;

               /*     ..Interpolate and accumulate */
               if (ir < NetNr) {

                  /*  ..Interpolate between adjacent samples */
                  alpha = ir - normRprime;

                  /*tmp = radpsf(ir)*alpha + radpsf(ir+1)*(1.0-alpha)*/
                  tmp = *(radpsf+ir-1) * ((float)ir-normRprime)
                        + *(radpsf+ir+1-1) * (normRprime+1.0-(float)ir);
                  }
               else {
                  /*This should disapear when we become confident
                    that the conservative guess for Maxir is OK */
                  if (ir > Maxir) {
                     Maxir = ir;
                     fprintf(stderr,"Found an ir > Maxir. \n");
                     fprintf(stderr,"Maxir set to %d\n",ir);
                     fprintf(stderr,"check the code that guesses Maxir again.\n");
	             appod = 1.0/(Maxir-Nr);
                     }

                  /*..Apodize from last available sample to zero 
            		to avoid sharp edges in the PSF */
                  /* ..Apodize by linear interpolation to zero at ir = Maxir */
                  weight = (Maxir-(float)ir)*appod;

                  tmp = *(radpsf+NetNr-1) * weight;
                  }  

               illum = illum + tmp;
            }   /*           end loop for n1 */
         }   /*           end loop for n2 */

         /* IMPORTANT
            The dimension of the array psfcond is now Mxy by Mxy,
            instead of (Pxy+2) by Pxy. */

         /* psfcond(ix,iy) = illum */
         *(psfcond+(ix-1)+(iy-1)*Mxy) = illum;

      } /*         end loop for ix */
   } /*         end loop for iy */

   /* ..Even symmetric replication into the other three quadrants */

   /* IMPORTANT
      The dimension of the array psfcond is now Mxy by Mxy,
      instead of (Pxy+2) by Pxy. 
   */

   evenrepr(psfcond, Mxy, Mxy);

   /*  ..Multiply illumination and objective psf (into psfobj) */
   mult3drm(psfobj, psfcond, Mxy*Mxy);
     
   /* sum neighboring pixels to downsample the BIG PSF */
   /* note dimension goes from Mxy to Lxy = Nxy*bin */
   Sum4N4NToNN(psfobj,psfbin1,Lxy,osamp); /* routine is in rotsum.c */

   /* Bin here */
   /* If binning by 1 copy psf to psfbin  */
   if (bin == 1) /* Lxy = Nxy */
     cp3dr (psfbin, psfbin1, Nxy*Nxy);
   else {
      /* Add pixels together */
      /* Clean psfbin array to use as accumulator */
      init3dr(psfbin, 0.0, Nxy*Nxy);

      /* binning  */
      for(ix=1;ix<=Lxy;ix++){
         jx = ((ix-1)/bin) + 1;
         for(iy=1;iy<=Lxy;iy++){
            jy = ((iy-1)/bin) + 1;
	    *(psfbin+(jx-1)+(jy-1)*Nxy) += *(psfbin1+(ix-1)+(iy-1)*Lxy);
            }
         }
      }

   WritePlane(iz-1,outpsf,psfbin,Nxy,Nxy,Nxy,Nxy);
   if(symmetric && (iz > 1))
       WritePlane(Nz-iz+1,outpsf,psfbin,Nxy,Nxy,Nxy,Nxy);


}  /* end of "for(iz=1;iz<=Zup;iz++)" */


keeplog (logfp, "DONE");
if(logfp != stderr) fclose(logfp);

printf("don't free up the memory\n");
/*free(radpsf);
free(psfbin);
free(psfobj);
free(psfcond);
free(pupil);*/

}



