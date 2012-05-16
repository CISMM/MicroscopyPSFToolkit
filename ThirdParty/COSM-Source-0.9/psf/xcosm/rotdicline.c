/**************************************************************

rotdicline.c

       Purpose:
            Creates a 3-D psf volume from a given XZ slice.
            The XZ slice MUST have been created by the program
            'gibson_xzslice' because of values inserted into the header.

            The slice was oversampled in gibson_xzslice, so this program
            reads in a row of that slice, radially sweeps it 360 degrees
            (centered at 0,0), takes DFT of the resulting 2d plane, then
            extracts the 4 corners in Fourier Space (corresponding to the
            removal of frequences > 1/(2*Deltax)),
            forces H(N/2+1,N-eta) = Conjugate{ H(N/2+1,eta) } 
            for eta = [1,N/2], changes the value at N/2,N/2 to 'real'
            by RE(i) = |i|*sign(RE(i)), IM(i) = 0.0,
            then takes reverse DFT and stores the result as a plane
            of the output psf. 

            The program also checks for a symmetric psf, in which case
            the above procedure need only be performed for N/2+1 rows        
            and then duplicated.

            The resulting psf is centered at 0,0,0.  

**********************************************************/

#include "misc.h"
#ifndef WIN32
#include <sys/param.h>
#endif
#include <math.h>
#include "washu.h"
 

extern char lognm[];

char tmpstring[255];

extern void fftn(float *r, float *i, unsigned nn[], int ndim, int isign);

#define PI      3.1415926535
#define DatOf1D(dat,x)                  *(dat + (x))
#define DatOf2D(dat,x,y,wid)            *(dat + (x) + (y)*wid)



int N,N2,N4,Nover2,Nz,Nr;

void MakeSumToOne(float *dat, int dnx, int dny, int dnz)
{
float dval,*fptr = dat;
double dsum = 0.0;
int i,tot_size = dnx*dny*dnz;
for(i=0;i<tot_size;i++)
 { dsum += (double)*fptr; fptr++; }
fptr = dat;
for(i=0;i<tot_size;i++)
 { dval = (double)*fptr; *fptr++ = (float)(dval/dsum); }
}

/* Convert complex number to real by Re(i) = |i|*sign(Re(i)), Im(i) = 0.0 */

void ChangeToReal(float *rpart, float *ipart, int x, int y, int wid)
{
   float rval,ival,mag;
   rval = DatOf2D(rpart,x,y,wid);
   ival = DatOf2D(ipart,x,y,wid);
   mag = sqrtf(rval*rval + ival*ival);
   if(rval < 0.0) mag = -mag;
   DatOf2D(rpart,x,y,wid) = mag;
   DatOf2D(ipart,x,y,wid) = 0.0;
}

/* Downsize 4N x 4N plane to N x N by taking four corners */
void Extract4N4NToNNDisplay(float *src)
{
   int x,y;
   float temp;
   int secondy,secondx;
   float *dptr;
   secondy = secondx = N4-Nover2+1;
   for(y=0;y<=Nover2;y++)
   {
      for(x=0;x<=Nover2;x++)
      {
        temp=DatOf2D(src,x,y,N4);
        if (temp)
        printf("Display DatOf2D:%lf",temp);
      }
      for(x=secondx;x<N4;x++)
      {
        
        temp=DatOf2D(src,x,y,N4);
        if (temp)
           printf("Display DatOf2D:%lf",temp);
      }
   }

   for (y=secondy;y<N4;y++)
   {
     for(x=0;x<=Nover2;x++)
     {
        temp=DatOf2D(src,x,y,N4);
        if (temp)
        printf("Display DatOf2D:%lf",temp);
     }
     for(x=secondx;x<N4;x++)
     {
        temp=DatOf2D(src,x,y,N4);
        if (temp)
        printf("Display DatOf2D:%lf",temp);
     }
   }
}

/* Downsize 4N x 4N plane to N x N by taking four corners */
void Extract4N4NToNN(float *src, float *dest)
{
   int x,y;
   int secondy,secondx;
   float *dptr;
   secondy = secondx = N4-Nover2+1;
   dptr = dest;
   for(y=0;y<=Nover2;y++)
   {
      for(x=0;x<=Nover2;x++)
        *dptr++ = DatOf2D(src,x,y,N4);
      for(x=secondx;x<N4;x++)
        *dptr++ = DatOf2D(src,x,y,N4);
   }

   for (y=secondy;y<N4;y++)
   {
     for(x=0;x<=Nover2;x++)
        *dptr++ = DatOf2D(src,x,y,N4);
     for(x=secondx;x<N4;x++)
        *dptr++ = DatOf2D(src,x,y,N4);
   }
}
void WritePlane2(int plane, float *dout, float *din, int w, int h)
{
     int i;
     int pls;
     float *outloc,*inloc,fval;
     pls = w*h;
     inloc = din;
     fflush(stdout);
     outloc = dout + plane*pls;
     //printf("WritePlane2 Nz-iz is %d,w is %d,h is %d,pls is %d,plane*pls is %d\n",plane,w,h,pls,(plane*pls));

     for (i=0;i<pls;i++)
     {
        fval = *inloc++;
        *outloc++ = fval;
     }
}

void rotdicline(float *outpsf_re, 
                float *outpsf_im,
                osm_ds *head_re, 
				osm_ds *head_im,
                int bin,
                float DeltaX, 
                float DeltaPhi, 
                float AmpRatio,
                int tandem)
{
   int num_pix;
   float twoPInum_pix;
   float sinfactor, cosfactor, ratfactor, angle;
   float f, deltaf, delta, fmax, fmax_half;

   FILE *fplog;
   int i,j,xi,eta,index;
   int symmetric;
   int oversmp;
   int ix,iy,ir,iz,ntot;
   int data_nx,data_ny,data_pl;
   int Nzwanted;
   unsigned nn[2];
  
   float fval,deltaxy,deltaxy1,deltar,r_hat,r,x,y;
   float *curptr,*rowdat,*outimg,*outptr,*fptr,*outimgosamp,*outptrosamp;
   float *outimg2,*outptr2;
   float *rowdat2;
   float *ptrr,*ptri;
   float *real4N4N,*imag4N4N;
   float *realNN,*imagNN;
   float interp_val,radval1,radval2,perc1,perc2;
   float interp_val2, temp;

   ratfactor = 1.0 - 2.0 * AmpRatio;

   rowdat = (float*)head_re->data;
   rowdat2 = (float*)head_im->data;
   data_nx = (int)(head_im->nx);
   data_ny = (int)(head_im->ny);
   data_pl = data_nx * data_ny;
   Nz = data_ny;
   Nr = data_nx;

   /* This is tentative, make sure gibson_xzslice places values in same place */
   oversmp = (int)(head_im->xstart);          /* oversmp */
   N = (int)(head_im->ystart);                /* Nxy */
   symmetric = (int)(head_im->zstart);        /* symmetr */
   ntot = N*N;
   N2 = N*2; N4 = N*4; Nover2 = N/2;
   /* Note: deltaxy and deltar are based on Nxy/2 = Nover2 */
   deltar = head_im->xlength;
   deltaxy = deltar*(float)oversmp;
   fmax = 1/ deltaxy ;
   fmax_half = fmax / 2.;
   deltaf = fmax / N;
   delta = 1 / N;
   fprintf(stderr,"DeltaX= %f , DeltaPHI= %f ,AmpRatio %f\n", DeltaX,DeltaPhi,AmpRatio);
   fprintf(stderr,"deltaxy= %f (mm), deltaf= %f (1/mm)\n", deltaxy,deltaf);

   if((fplog=fopen(lognm,"w"))==(FILE*)NULL)
   {
      fprintf(stderr,"ERROR: Can't open '%s' for write/append\n",lognm);
      exit(1);
   }

   fprintf(fplog,"Shear along x in image space (mm) : %f\n",DeltaX);
   fprintf(fplog,"                Phase bias (rads) : %f\n",DeltaPhi);
   fprintf(fplog,"                  Amplitude ratio : %f\n",AmpRatio);
   fprintf(fplog,"   Pixel size in image space (mm) : %f\n",deltaxy);
   fprintf(fplog," Pixel size in freq. space (1/mm) : %f\n",deltaf);
   fprintf(fplog,"\n");
   fclose(fplog);

   DeltaPhi = DeltaPhi / 2.0;
   DeltaX = DeltaX / 2.0;
   num_pix = (int)(DeltaX / deltaxy + 0.5);
   twoPInum_pix = 2.*PI*num_pix;

   fprintf(stderr,"DeltaX/2= %f mm (num_pix=%d) DeltaPhi/2=%f rads\n",DeltaX, num_pix ,DeltaPhi);

   /* New deltaxy for 4x oversampling */
   deltaxy1 = deltaxy/4.0;

   if (symmetric == 0)
      Nzwanted = Nz;
   else
      Nzwanted = Nz/2+1;

   /* Output psf data */
   if ((outimg=(float*)calloc(sizeof(float),N*N*Nz))==(float*)NULL)
   {
      fprintf(stderr,"ERROR: Out of memory (1).\n");
      exit(1);
   }
   if ((outimg2=(float*)calloc(sizeof(float),N*N*Nz))==(float*)NULL)
   {
     fprintf(stderr,"ERROR: Out of memory (1).\n");
     exit(1);
   }

   /* Real part of the 4Nx4N oversampled 2-D plane */
   if((real4N4N=(float*)calloc(sizeof(float),N4*N4))==(float*)NULL)
   {
     fprintf(stderr,"ERROR: Out of memory (2).\n");
     exit(1);
   }
   if ((imag4N4N=(float*)calloc(sizeof(float),N4*N4))==(float *)NULL)
   {
      fprintf(stderr,"ERROR: Out of memory (3).\n");
      exit(1);
   }

   /* Real part of the NxN 2-D plane */
   if ((realNN=(float*)calloc(sizeof(float),N*N))==(float *)NULL)
   {
     fprintf(stderr,"ERROR: Out of memory (6).\n");
     exit(1);
   }

   /* Imaginary part of the NxN 2-D plane */
   if ((imagNN=(float*)calloc(sizeof(float),N*N))==(float *)NULL)
   {
     fprintf(stderr,"ERROR: Out of memory (7).\n");
     exit(1);
   }

   outptr = outimg;
   outptr2 = outimg2;

   for(iz=0;iz<Nzwanted;iz++)
 {
 /* Sweep psf radial data into 1 plane */
 for(iy=0;iy<=N2;iy++)
  {
  y = (float)iy * deltaxy1;
  for(ix=0;ix<=N2;ix++)
   {
   x = (float)ix * deltaxy1;
   r = sqrtf(x*x + y*y);
   r_hat = r/deltar;
   ir = (int)r_hat;
   radval1 = DatOf1D(rowdat,ir);
   radval2 = DatOf1D(rowdat,ir+1);
   perc2 = r_hat - (float)ir;
   perc1 = 1.0 - perc2;
   interp_val = radval1*perc1+radval2*perc2;

   radval1 = DatOf1D(rowdat2,ir);
   radval2 = DatOf1D(rowdat2,ir+1);
 interp_val2 = radval1*perc1+radval2*perc2;
   //printf("interp_val is %lf\n",interp_val);
   //printf("interp_val2 is %lf\n",interp_val2);

   /* Upper left */
   if((iy < N2)&&(ix < N2))
        {
     DatOf2D(real4N4N,ix,iy,N4) = interp_val;
     DatOf2D(imag4N4N,ix,iy,N4) = interp_val2;
        }
   /* Upper right */
   if(ix > 0)
        {
     DatOf2D(real4N4N,N4-ix,iy,N4) = interp_val;
     DatOf2D(imag4N4N,N4-ix,iy,N4) = interp_val2;
        }
   /* Lower left */
   if(iy > 0)
        {
     DatOf2D(real4N4N,ix,N4-iy,N4) = interp_val;
     DatOf2D(imag4N4N,ix,N4-iy,N4) = interp_val2;
        }
   /* Lower right */
   if((ix > 0)&&(iy > 0))
        {
     DatOf2D(real4N4N,N4-ix,N4-iy,N4) = interp_val;
     DatOf2D(imag4N4N,N4-ix,N4-iy,N4) = interp_val2;
        }
   }
  }



 /* Now take 2D forward DFT of 4Nx4N array calculated above */
/* forward here means using exp(-j...) and that's why I use -1 in the call
which is really the reverse for the definition of this fft routine */
 nn[0] = N4; nn[1] = N4;
 fftn(real4N4N,imag4N4N,nn,2,-1);

 /*  
   Now, downsample in X/Y direction by 1/4th
   removing unwanted frequencies > 1/(2deltax) 
 */

 //Extract4N4NToNNDisplay(real4N4N);

 Extract4N4NToNN(real4N4N,realNN);
 Extract4N4NToNN(imag4N4N,imagNN);

/* ChangeToReal(realNN,imagNN,Nover2,Nover2,N); */
printf("twoPInum_pix=%f, DeltaPhi=%f",twoPInum_pix,DeltaPhi);
if (ratfactor == 0.0) {
  for(i=0;i<N;i++)
  {
    ir= i * N;
    for (j=0;j<N;j++)
    {
       index=ir + j;
       f = (float) j * deltaf; 
       if (f > fmax_half)   
          f -= fmax;
       f = f/fmax;

       /* multiply by phase factor for DIC optics */
       /* phase factor was changed on 1/31/95; also on 7/30/96 I added the "-"
          do be consistent with the FFT exp sign. */
       /* on 11/23/96 I realized that the "-" here was a mistake because I call
          the fft with the opposite sign */
       sinfactor =  sin(twoPInum_pix * f + DeltaPhi) ;
       //printf("sinfactor is %lf\n",sinfactor);
       temp =  sinfactor * imagNN[index];
       imagNN[index] = - sinfactor * realNN[index];
       realNN[index] = temp;
     }
   }
 } 
 else {
 	for(i=0;i<N;i++)
	  {
            ir= i * N;
  	    for (j=0;j<N;j++)
   	    {
   	       index=ir + j;
   	       f = (float) j * deltaf; 
   	       if (f > fmax_half)	
                  f -= fmax;
	       f = f/ fmax;
   	       /* multiply by phase factor for DIC optics */
   	          angle = twoPInum_pix * f + DeltaPhi;
	       /* note that the rest has been changed as of 3/8/95 */
   	       sinfactor = sin(angle);
   	       cosfactor = cos(angle) * ratfactor;
   	       temp = cosfactor * realNN[index] + sinfactor * imagNN[index];
   	       imagNN[index] = cosfactor*imagNN[index] - sinfactor*realNN[index] ;
   	       realNN[index] = temp;
    	    }
         }
    }

 /* Take inverse transform to get plane of the PSF */
 nn[0] = N; nn[1] = N;
 fftn(realNN,imagNN,nn,2,1);

/*Scale by a factor of 4x4=16 to accomodate for the FFT scale factor and the
   downsampling by a factor of 4 in x and y */
 ptrr = realNN; ptri = imagNN;
 for(ir=0;ir<ntot;ir++)
	 {
	  *ptrr++ /= 16.0;
	  *ptri++ /= 16.0;
	  }

 /* Copy to 3-D PSF */

 printf("\nReal creation\n");
 WritePlane2(iz,outpsf_re,realNN,N,N);
 if((symmetric != 0)&&(iz > 0)&&(iz < Nzwanted-1))
	WritePlane2(Nz-iz,outpsf_re,realNN,N,N);

 printf("Imaginary creation\n");
 WritePlane2(iz,outpsf_im,imagNN,N,N);
 if((symmetric != 0)&&(iz > 0)&&(iz < Nzwanted-1))
	WritePlane2(Nz-iz,outpsf_im,imagNN,N,N);

 /* Increment row counter to point to next row */
 rowdat += Nr;
 rowdat2 += Nr;
 }
 
fprintf(stderr,"DONE.\n");


/* Not needed, as gen_eig normalizes OTF for you 
MakeSumToOne(outimg,N,N,Nz);

*/
/*
head->data = (char*)outimg;
head->nx = (headtype)N;
head->ny = (headtype)N;
head->nz = (headtype)Nz;
head->user15 = TVAL;
head->user16_footer_size = 0;
//update_stats(head);

sprintf(tmpstring,"%s.re",argv[2]);
write_dataset(head,tmpstring);


head->data = (char*)outimg2;
update_stats(head);
sprintf(tmpstring,"%s.im",argv[2]);
write_dataset(head,tmpstring);

*/ 


}

