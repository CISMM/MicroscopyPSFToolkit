/***************************************************************************
  COPYRIGHT 1996 BIOMEDICAL COMPUTER LABORATORY, WASHINGTON UNIVERSITY, MO
****************************************************************************/
#ifndef _FFT3D_H
#define _FFT3D_H
typedef struct {
	float re;
	float im;
	} fcomplex;

extern void real_fft3d(float *cplxdat, int nx, int ny, int nz, int xdim, int ydim, int dir);

#ifdef SGI_FFT

#define FORWARD -1
#define REVERSE  1
#define FACTOR_SPACE	15

extern void sscal3d( int nx, int ny, int nz, float alpha, float *y, int ld1,int ld2);
extern float *sfft3dui( int n1, int n2, int n3, float *save);
extern int sfft3du( int job, int n1, int n2, int n3, float *array, int ld1, int ld2, float *save);

#else
#define FORWARD  1
#define REVERSE -1
#endif

#endif
