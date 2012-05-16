/***************************************************************************
  COPYRIGHT 1996 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "fft3d.h"


#ifdef _MISCC_
#define EXTERN 
#else
#define EXTERN extern
#endif

EXTERN float amax1(float a, float b);
EXTERN float amin1(float a, float b);
EXTERN float myfabs(float f);
EXTERN void init3dr(float *array, float mconst, int siz);
EXTERN void cp3dr(float *dest, float *orig, int siz);
EXTERN void idefaulter(char *prompt, int *var);
EXTERN void defaulter(char *prompt, float *var);
EXTERN void getpupil(float *pupil, float diam, int nx, int ny);
EXTERN void evenrepr(float *array, int nx, int ny);
EXTERN int  intlog2(int number);
EXTERN void keeplog(FILE *fp, char *message);
EXTERN void add2columns(float *pt, int nx, int nyz);
EXTERN void mult3dcm(fcomplex *array1, fcomplex *array2, int siz);
EXTERN void mult3drm(float *array1, float *array2, int siz);
EXTERN void norm3dcm(fcomplex *array, float *peak, int siz);
EXTERN void norm3drm(float *array, float *peak, int siz);
EXTERN void r2xy(float *fxy, float *fr, int nx, int ny, int nr, float ratio);

