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

#include <blitz/array.h>
#include <complex>
#include <fftw3.h>
#include <ostream.h>
#include "fftwInterface.h"

using namespace blitz;
using namespace cosm;

int main(
){
     int i;
     const int N = 16;
     int size = 2;
     int n[1]  = {N};

     // first fftw api directly using fftw data types
     cout << "fftw api using fftw data types"<<endl;
     double* in;
     fftw_complex *out;
     fftw_plan p;
     in = (double*) fftw_malloc(sizeof(double)*N);
     out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*((N/2)+1));
     // forward fft
     p = fftw_plan_dft_r2c(1, n, in, out, FFTW_ESTIMATE);
     cout << "in: ";
     for ( i = 0; i < size; i++ ) {
	in[i] = 1;
        cout << in[i] <<" ";
     }
     for ( i = size; i < N; i++ ) {
	in[i] = 0.0;
        cout << in[i] <<" ";
     }
     cout << endl;
     fftw_execute(p);
     cout << "out: ";
     for ( i = 0; i <  N/2+1; i++ ) {
        cout <<"("<< (out[i])[0]<<","<<(out[i])[1]<<") ";
     }
     cout << endl;
     fftw_destroy_plan(p);
     // reverse fft
     p = fftw_plan_dft_c2r(1, n, out,in, FFTW_ESTIMATE);
     fftw_execute(p);
     cout << "in(FFTed): ";
     for ( i = 0; i < N; i++ ) {
        cout << in[i]/N <<" ";
     }
     cout << endl;
     fftw_destroy_plan(p);
     fftw_free(in);
     fftw_free(out);

     // use fftw api directly using Array data types
     cout << "fftw api using Array data types"<<endl;
     Array<double,1> A(N);
     Array<std::complex<double>,1> B(N/2+1);
     Array<double,1> C(N);
     firstIndex a;
     A = where(a < size, 1, 0);
     cout <<"A: "<<A<<endl;
     // forward fft
     p = fftw_plan_dft_r2c(1, n, A.data(), reinterpret_cast<fftw_complex*>(B.data()), FFTW_ESTIMATE);
     fftw_execute(p);
     cout <<"B: "<<B<<endl;
     fftw_destroy_plan(p);
     // reverse fft
     p = fftw_plan_dft_c2r(1, n, reinterpret_cast<fftw_complex*>(B.data()), C.data(), FFTW_ESTIMATE);
     fftw_execute(p);
     C /= N;
     cout <<"C: "<<C<<endl;

     // use fftwInterface using Array data types
     cout << "fftwInterface using Array data types"<<endl;
     fftwInterface<double,1> fftw;
     cout <<"A: "<<A<<endl;
     // forward fft
     fftw.plan(A,B);
     fftw.execute();
     cout <<"B: "<<B<<endl;
     // reverse fft
     fftw.plan(B,C);
     fftw.execute();
     C /= N;
     cout <<"C: "<<C<<endl;

     // use fftwInterface global functions 
     cout << "use fftwInterface global functions"<<endl;
     // forward fft
     cout <<"A: "<<A<<endl;
     Array<std::complex<double>,1> D = forwardFFT(A);
     cout <<"B: "<<B<<endl;
     // reverse fft
     Array<double,1> E = inverseFFT(D);
     cout <<"C: "<<C<<endl;
      
     return 0;
}
