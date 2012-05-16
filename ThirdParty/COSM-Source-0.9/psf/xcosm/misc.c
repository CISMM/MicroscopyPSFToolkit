/***************************************************************************
  COPYRIGHT 1996 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
#include "misc.h"

/***************************************************************
PURPOSE:
To initialize a 3-D real array to a constant

PARAMETERS:
array  (float, OUTPUT) destination array 
constant  (float, INPUT), Constant to assign to all element of array
Ntot    size of array
***************************************************************/

void init3dr(float *array, float fconstant, int Ntot)
{
int	i;
float	*tmp;
tmp=array;
for(i=0;i<Ntot;i++) *tmp++ = fconstant;
}


float myfabs(float f)
{
if(f < 0.0) return(-f);
else return(f);
}

float amax1(float a, float b)
{
if(a > b) return(a);
else      return(b);
}

float amin1(float a, float b)
{
if(a > b) return(b);
else      return(a);
}

/*************************************************************
PURPOSE:
To copy 3-D real array origin into 3-D real array dest.

PARAMETERS:
dest  (float, OUTPUT) destination array
origin: (float, INPUT) original array 
Ntot  size of array, =Nx*Ny*Nz
***************************************************************/

void cp3dr(float *dest, float *origin, int Ntot)
{
int	i;
float	*tmp1,*tmp2;

tmp1=dest;
tmp2=origin;
for(i=0;i<Ntot;i++) *tmp1++ = *tmp2++;
}


/*********************************************************
PURPOSE:
	Prompts for the value of an integer variable
	The variable must take only non-negative values

	User is promped for a quantity, and a default is printed
	together with the prompt string.  If the user types in zero
	or a negative number the default value is taken, other entries
	are assigned to 'variable'
PARAMETERS:
	prompt: (INPUT) string to be used to prompt
		for the quantity.
	variable (int;INPUT/OUTPUT) On input contains the default
		value of the variable.  If the variable can be zero the
		default MUST be zero, otherwise there is no way to 
		make it take that value.  On output it has the value
		given by the user.
************************************************************/
void idefaulter (char *prompt, int *variable)
{
int	tmp;

fprintf(stderr,"%s\n (default=%8d):\n",prompt, *variable);
fscanf(stdin,"%d",&tmp);
if (tmp > 0) *variable = tmp;
else if(tmp <0) 
   fprintf(stderr,"your input is negative; default is used instead.\n");
}


/*********************************************************
PURPOSE:
	Prompts for the value of a float variable
	The variable must take only non-negative values

	User is promped for a quantity, and a default is printed
	together with the prompt string.  If the user types in zero
	or a negative number the default value is taken, other entries
	are assigned to 'variable'
PARAMETERS:
	prompt: (INPUT) string to be used to prompt
		for the quantity.
	variable (float;INPUT/OUTPUT) On input contains the default
		value of the variable.  If the variable can be zero the
		default MUST be zero, otherwise there is no way to 
		make it take that value.  On output it has the value
		given by the user.
************************************************************/
void defaulter (char *prompt, float *variable)
{
float	tmp;

fprintf(stderr,"%s\n (default=%f):\n",prompt, *variable);
fscanf(stdin,"%f",&tmp);
if (tmp > 0.0) *variable = tmp;
else if(tmp <0.0) 
   fprintf(stderr,"your input is negative; default is used instead.\n");
}


/**********************************************************
PURPOSE:
to sample a circular pupil, replicate with circular symmetry and
PARAMETERS:
pupil (OUPUT, flota 2-D ARRAY) array that will contain the Fourier
transform of the pupil.  This subroutine uses it as an array
of reals.
diam: (INPUT, int) Pupil diameter in pixels
Nnx, Nny: (INPUT, int) dimensions of array pupil
Array pupil is treated like a 1-D array. But only a part of
the array is explicitly calculated in this function. the rest of the
elements are obtained by symmetry argument. Since Nnx and
Nny can be different, we need two indices in this function,
instead of only one in other cases.
If all elements are explicitly calculated in this function,
we would need only one index Ntot.
***********************************************************/

void getpupil (float *pupil, float diam, int Nnx, int Nny)
{
int	ix, iy, HalfNnx, HalfNny;
/*   ..Square of pupil radius  */
float	radsq;
/*   ..running coordinates */
float	x, y, ysq, rd;


HalfNny = Nny/2;
HalfNnx = Nnx/2;
radsq = diam*diam/4.0;

/*   ..Clearing pupil array  */
init3dr(pupil, 0.0, Nnx*Nny);
for(iy=1;iy<=HalfNny+1;iy++){
   y = (float) iy-1;
   ysq = y*y;
   for(ix=iy;ix<=HalfNnx+1;ix++){
      x = (float) ix-1;
      rd = x*x + ysq;
      if (rd < radsq) {
         *(pupil+(ix-1)+(iy-1)*Nnx) = 1.0;
         *(pupil+(iy-1)+(ix-1)*Nnx) = 1.0;
      }
   }
}

/*   ..Replicate to the other three quadrants */
evenrepr (pupil, Nnx, Nny);
}


/********************************************************
PRPOSE:
To do a periodic replication of a portion of a real, 2-D array from
ix=1, Nnx/2, iy = 1, Nny/2
into the rest of an Nnx by Nny array assuming even symmetry in 
both ix and iy.

PARAMETERS:
array (INPUT/OUTPUT, 2-D float ARRAY) 
Nnx, Nny: dimensions of the array, x-index runs faster than y-index
********************************************************/
void evenrepr(float *array, int Nnx, int Nny)
{
/* __Internal variables */
int	HalfNnx, HalfNny;
int	ix, iy;

HalfNnx = Nnx/2;
HalfNny = Nny/2;
/*   ..Replicate quadrant I into II */
for(iy=1;iy<=HalfNny+1;iy++){
   for(ix=2;ix<=HalfNnx;ix++){
/*    array(Nnx-ix+2, iy) = array(ix, iy)   */
      *(array+(Nnx-ix+2-1)+(iy-1)*Nnx) = *(array+(ix-1)+(iy-1)*Nnx);
   }
}
/*   ..Replicate quadrant I into IV and II into III */
for(ix=1;ix<=Nnx;ix++){
   for(iy=2;iy<=HalfNny;iy++)
/*    array(ix, Nny-iy+2) = array (ix,iy) */
      *(array+(ix-1)+(Nny-iy+2-1)*Nnx) = *(array+(ix-1)+(iy-1)*Nnx);
}
}


/******************************************************
Return the largest integer less than or equal to the
base-two logarithm ofthe integer 'number'.
*******************************************************/
int intlog2(int number)
{
double ln2=0.69314718056; /* natural log of 2 */
return (int) (log((double) number + 0.5) / ln2);
}


/*************************************************
Append 'message' to end of logfile already opened in
calling program. The message is written together
with the local date and time. The logfile is NOT
closed upon exiting 'keeplog'.
**************************************************/
void keeplog(FILE *fp, char *message)
{
time_t	t;

t=time(NULL);
fprintf(fp,"%s\t%s\n",asctime(localtime(&t)),message);
fflush(fp);
}


/********************************************************
insert two float spaces (not used) after every Nnx elements
in float array pt. The effect is to strech pt by 2*Nnyz
elements. Array pt must have enough space for these
extra elements.
********************************************************/
void add2columns(float *pt, int Nnx, int Nnyz)
{
int	ix,iyz;
float	*tmp1,*tmp2;

for(iyz=Nnyz;iyz>0;iyz--){
   tmp1=pt+(Nnx+2)*iyz-3;
   tmp2=pt+Nnx*iyz-1;
   for(ix=Nnx;ix>0;ix--) *tmp1-- = *tmp2--;
}
}



/********************************************************
remove two float spaces (not used) after every Nnx elements
in float array pt. The effect is to squeeze pt by 2*Nnyz
elements. 
********************************************************/
void remove2columns(float *pt, int Nnx, int Nnyz)
{
int	ix,iyz;
float	*tmp1,*tmp2;

for(iyz=0;iyz<Nnyz;iyz++){
   tmp1=pt+Nnx*iyz;
   tmp2=pt+(Nnx+2)*iyz;
   for(ix=0;ix<Nnx;ix++) *tmp1++ = *tmp2++;
}
}



/**************************************************************
PURPOSE 
	multiply two 3-D arrays of complex numbers element by element.
	the result is stored in the first array, i.e. 
	array1(ix,iy,iz) = array1(ix,iy,iz)*array2(ix,iy,iz) 
		for ix = 1, ... Nx; iy = 1, ... , Ny; iz = 1, ... Nz
        Ntot=Nx*Ny*Nz
        The storing order of the 3D arrays is determined in the 
	calling program. As far as mult3dcm is concerned, they're
	1D arrays with dimension Ntot.
PARAMETERS:
   array1: (fcomplex/INPUT, OUTPUT) pointer to array1;
	on output it contains the product.
   array2: (fcomplex/INPUT) pointer to the other 
        array in the product.  The array is unchanged on output.
   Ntot(int/INPUT) size of the two arrays
****************************************************************/
void mult3dcm(fcomplex *array1, fcomplex *array2, int Ntot)
{
int		i;
fcomplex	*pt1,*pt2;
float		tmp1,tmp2;
pt1=array1;
pt2=array2;
for(i=0;i<Ntot;i++){
   tmp1 = pt1->re * pt2->re - pt1->im * pt2->im;
   tmp2 = pt1->re * pt2->im + pt1->im * pt2->re;
   pt1->re = tmp1;
   pt1->im = tmp2;
   pt1++;
   pt2++;
}
}


/**************************************************************
PURPOSE 
	multiply two 3-D arrays of float numbers element by element.
	the result is stored in the first array, i.e. 
	array1(ix,iy,iz) = array1(ix,iy,iz)*array2(ix,iy,iz) 
		for ix = 1, ... Nx; iy = 1, ... , Ny; iz = 1, ... Nz
        Ntot=Nx*Ny*Nz
        The storing order of the 3D arrays is determined in the 
	calling program. As far as mult3drm is concerned, they're
	1D arrays with dimension Ntot.
PARAMETERS:
   array1: (float/INPUT, OUTPUT) pointer to array1;
	on output it contains the product.
   array2: (float/INPUT) pointer to the other 
        array in the product.  The array is unchanged on output.
   Ntot(int/INPUT) size of the two arrays
****************************************************************/
void mult3drm(float *array1, float *array2, int Ntot)
{
int	i;
float	tmp,*pt1,*pt2;
pt1=array1;
pt2=array2;
for(i=0;i<Ntot;i++){
   tmp = *pt1 * *pt2;
   *pt1 = tmp;
   pt1++;
   pt2++;
}
}





/************************************************************
To find the maximum modulus of a complex 3-D array and divide each 
element by this maximum.
Array is treated like a 1-D array in this function. The storing
order of array is determined by the calling program.
Array size is Ntot=N1*N2*N3. 
On output 'array' is normalized to maximum modulus = 1.
If maximum = 0, array is not divided.

PARAMETERS:
array:			(INPUT/OUTPUT) pointer to complex 3-dimensional
			array to normalize
			On output contain the normalized array
peak:	  		(INPUT/OUTPUT) if on input peak = 0.0
			the array is normalized by the maximum 
			absolute value 	found in the array and 
			'peak' set to this value,
			if on input peak .ne. 0 each element of 
			the array is divided by 'peak'. 
			In this case the variable peak doesn't change
Ntot:			(INPUT) array dimension 
*********************************************************************/
void norm3dcm ( fcomplex *array, float *peak, int Ntot )
{
int	i;
/*****    ..Criterium to regard somethig as almost zero  ****/
float	EPSILON = 1.0e-30;
float	tmp;
fcomplex	*pt;

pt=array;

if (*peak == 0.0)
   for(i=0;i<Ntot;i++){
      tmp = (float) sqrt(pt->re * pt->re + pt->im * pt->im); 
      pt++;
      *peak = (tmp > *peak) ? tmp : *peak;
   }

if (*peak < EPSILON) {
   fprintf(stderr,"ERROR IN ROUTINE NORM3DCM MAX VALUE TOO SMALL\n");
   fprintf(stderr,"PEAK=%e\n",*peak);
   *peak = 0.0;
   return;
}

pt=array;
for(i=0;i<Ntot;i++){
   pt->re /= *peak;
   pt->im /= *peak;
   pt++;
}

}



/************************************************************
To find the maximum modulus of a real 3-D array and divide each 
element by this maximum.
Array is treated like a 1-D array in this function. The storing
order of array is determined by the calling program.
Array size is Ntot=N1*N2*N3. 
On output 'array' is normalized to maximum modulus = 1.
If maximum = 0, array is not divided.

PARAMETERS:
array:			(INPUT/OUTPUT) pointers to real 3-dimensional 
			real array to normalize;
			On output contain the normalized array
peak:	  		(INPUT/OUTPUT) if on input peak = 0.0
			the array is normalized by the maximum 
			absolute value 	found in the array and 
			'peak' set to this value,
			if on input peak .ne. 0 each element of 
			the array is divided by 'peak'. 
			In this case the variable peak doesn't change
Ntot:			(INPUT) array dimension 
*********************************************************************/
void norm3drm ( float *array, float *peak, int Ntot )
{
int	i;
/*****    ..Criterium to regard somethig as almost zero  ****/
float	EPSILON = 1.0e-30;
float	*pt;

pt=array;

if (*peak == 0.0)
   for(i=0;i<Ntot;i++){
      *peak = (*pt > *peak) ? *pt : *peak;
      pt++;
   }

if (*peak < EPSILON) {
   fprintf(stderr,"ERROR IN ROUTINE NORM3DCM MAX VALUE TOO SMALL\n");
   fprintf(stderr,"PEAK=%e\n",*peak);
   *peak = 0.0;
   return;
}

pt=array;
for(i=0;i<Ntot;i++){
   *pt /= *peak;
   pt++;
}

}


/************************************************************
PURPOSE:
take a sampled radially symmetric function fr(r) and compute
the correspoinding function over xy by bylinar interpolation
PARAMETERS:
fxy: (OUTPUT, 2-D float ARRAY) will contain fxy(x,y)
the output function is defined in four quadrants
(with periodic replication of quadrants II - IV)
fr : (INPUT, 1-D float ARRAY) contains the sampled function
the input function is sampled over non-negative r
Nnx, Nny: (INPUT, int) dimensions of array fxy
Nnr : (INPUT, int) Number of points in fr(r) (dimension of array fr)
ratio: (INPUT, float)  ratio of sampling distances Dxy/Dr, where
Dxy is the desired sampling in fxy(x,y) and Dr the available
sampling in fr(r).  ratio must be strictly positive
************************************************************/
void r2xy(float *fxy, float *fr, int Nnx, int Nny, int Nnr, float ratio)
{
int	iy, ix, ir;
float	x, y, rd, ysq;
/*   ..A number between zero and 1.0 for the interpolation */
float	alpha;
int	HalfNnx, HalfNny;
      
HalfNnx = Nnx/2;
HalfNny = Nny/2;

#ifdef DEBUG
printf("in r2xy: Nnx= %d, Nny=%d, Nnr = %d, HalfNnx=%d, ratio=%f\n",Nnx,Nny,Nnr, HalfNnx, ratio);
#endif

init3dr (fxy, 0.0, Nnx*Nny);
/*   ..Interpolate in the first quadrant */
for(iy=1;iy<=HalfNny+1;iy++){
   y = (float) iy-1;
   ysq = y*y;
   for(ix=iy;ix<=HalfNnx+1;ix++){
      x = (float) ix-1;
      rd = ratio * (float) sqrt(x*x+ysq);
      ir = (int) rd+1;
      if (ir < Nnr) { 
         alpha = ir-rd;
/* some debug stuff, probably will never be used */
         if (alpha > 1.0) {
            fprintf(stderr,"Error: alpha larger than 1.0 in r2xy\n");
            fprintf(stderr,"(ix,iy,alpha)=%d,%d,%f\n", ix, iy, alpha);
            exit(1);
         }
         else if(alpha < 0.0) {
            fprintf(stderr,"Error: alpha smaller than 0.0 in r2xy\n");
            fprintf(stderr,"(ix,iy,alpha)=%d,%d,%f\n", ix, iy, alpha);
            exit(1);
         }
/*       fxy(ix,iy) = fr(ir)*alpha + fr(ir+1)*(1.0-alpha)
         fxy(iy,ix) = fxy(ix,iy) */
         *(fxy+(ix-1)+(iy-1)*Nnx) = *(fr+(ir-1)) * alpha
 				 + *(fr+(ir+1-1)) * (1.0-alpha);
         *(fxy+(iy-1)+(ix-1)*Nnx) = *(fxy+(ix-1)+(iy-1)*Nnx);
      }
      else {
         fprintf(stderr,"Error: ir larger than Nnr\n");
         fprintf(stderr,"(ix,iy,ir,Nnr=)%d,%d,%d,%d\n",ix,iy,ir,Nnr);
         exit(1);
      }
   }
}
/*   ..Replicate to the other three quadrants */
evenrepr(fxy, Nnx, Nny);
}


