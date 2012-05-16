/****************************************************************************
 * Copyright (c) 2007 Einir Valdimarsson and Chrysanthe Preza
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 ****************************************************************************/

#ifndef _FFTW_INTERFACE_H
#define _FFTW_INTERFACE_H

#include <blitz/array.h>
#include <complex>
#include <fftw3.h>

using namespace blitz;
namespace cosm {

// forward declarations
class kindType;

template <typename T, int N>
class fftwInterface {

  public:

     // constructor
     fftwInterface();

     // destructor
     ~fftwInterface();

     // create a complex to complex plan;
     int plan( 
	const Array<std::complex<T>, N> &in, 
	Array<std::complex<T>, N> &out, 
	int sign,
	unsigned flags = FFTW_ESTIMATE
     );

     // create a real to complex plan;
     int plan( 
	const Array<T, N> &in, 
	Array<std::complex<T>, N> &out, 
	unsigned flags = FFTW_ESTIMATE
     );

     // create a complex to real plan;
     int plan( 
	const Array<std::complex<T>, N> &in, 
	Array<T, N> &out, 
	unsigned flags = FFTW_ESTIMATE
     );

     // create a real to real plan;
     int plan( 
	const Array<T, N> &in, 
	Array<T, N> &out, 
	kindType* kind,
	unsigned flags = FFTW_ESTIMATE
     );

     // execute the plan
     int execute();

  private:

     // wrapper support for different types

    fftwf_plan planC2C(
	int rank, int* n, complex<float>* in, complex<float>* out, int sign, unsigned flags
    ) { 
	return fftwf_plan_dft( rank, n, reinterpret_cast<fftwf_complex*>(in), reinterpret_cast<fftwf_complex*>(out), sign, flags); };
    fftw_plan planC2C(
	int rank, int* n, complex<double>* in, complex<double>* out, int sign, unsigned flags
    ) { 
	return fftw_plan_dft(rank, n, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), sign, flags); };
    fftwl_plan planC2C(
	int rank, int* n, complex<long double>* in, complex<long double>* out, int sign, unsigned flags
    ) { 
	return fftwl_plan_dft(rank, n, reinterpret_cast<fftwl_complex*>(in), reinterpret_cast<fftwl_complex*>(out), sign, flags); };


    fftwf_plan planR2C( 
	int rank, int* n, float* in, complex<float>* out, unsigned flags
    ) { 
	return fftwf_plan_dft_r2c(rank, n, in, reinterpret_cast<fftwf_complex*>(out), flags); };
    fftw_plan planR2C( 
	int rank, int* n, double* in, complex<double>* out, unsigned flags
    ) { 
	return fftw_plan_dft_r2c(rank, n, in, reinterpret_cast<fftw_complex*>(out), flags); };
    fftwl_plan planR2C( 
	int rank, int* n, long double* in, complex<long double>* out, unsigned flags
    ) { 
	return fftwl_plan_dft_r2c(rank, n, in, reinterpret_cast<fftwl_complex*>(out), flags); };

    fftwf_plan planC2R(
        int rank, int *n, complex<float>* in, float* out, unsigned flags
    ) { 
	return fftwf_plan_dft_c2r(rank, n, reinterpret_cast<fftwf_complex*>(in), out, flags); };
    fftw_plan planC2R(
        int rank, int *n, complex<double>* in, double* out, unsigned flags
    ) { 
	return fftw_plan_dft_c2r(rank, n, reinterpret_cast<fftw_complex*>(in), out, flags); };
    fftwl_plan planC2R(
        int rank, int *n, complex<long double>* in, long double* out, unsigned flags
    ) { 
	return fftwl_plan_dft_c2r(rank, n, reinterpret_cast<fftwl_complex*>(in), out, flags); };


    fftwf_plan planR2R(
        int rank, int* n, float* in, float* out, int* kind, unsigned flags
    ) { return NULL; }; 	// not implemented 
    fftw_plan planR2R(
        int rank, int* n, double* in, double* out, int* kind, unsigned flags
    ) { return NULL; }; 	// not implemented 
    fftwl_plan planR2R(
        int rank, int* n, long double* in, long double* out, int* kind, unsigned flags
    ) { return NULL; }; 	// not implemented 

    void destroyPlan( float tmp ) { 
	fftwf_destroy_plan((fftwf_plan)fftwPlan); 
    };
    void destroyPlan( double tmp ) { 
	fftw_destroy_plan((fftw_plan)fftwPlan); 
    };
    void destroyPlan( long double tmp ) { 
	fftwl_destroy_plan((fftwl_plan)fftwPlan); 
    };

    void executePlan( float tmp ) { 
	fftwf_execute((fftwf_plan)fftwPlan); 
    };
    void executePlan( double tmp ) { 
	fftw_execute((fftw_plan)fftwPlan); 
    };
    void executePlan( long double tmp ) { 
	fftwl_execute((fftwl_plan)fftwPlan); 
    };

  private:

    // not implemented
    fftwInterface( 
	const fftwInterface& 
    );
    fftwInterface& operator=( 
	const fftwInterface&  
    );

  private:
  
    void* fftwPlan;

};

template<typename T, int N>
Array<std::complex<T>,N> forwardFFT(
     const Array<T,N>& A
){
     fftwInterface<T,N> fftw;
     TinyVector<T, N> shape(A.shape()); 
     shape(N-1) = shape(N-1)/2+1;
     Array<std::complex<T>,N> B(shape);
     fftw.plan(A,B);
     fftw.execute();
     return B; 
};

template<typename T, int N>
Array<T,N> inverseFFT(
     const Array<std::complex<T>,N>& A
){
     fftwInterface<T,N> fftw;
     TinyVector<T, N> shape(A.shape()); 
     shape(N-1) = (shape(N-1)-1)*2;
     Array<T,N> B(shape);
     fftw.plan(A,B);
     fftw.execute();
     B /= B.size();
     return B; 
};

}

#include "blitz/fftwInterface.c"

#endif
