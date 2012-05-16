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
#include <fftw3.h>

namespace cosm {

// constructor
template<typename T, int N>
fftwInterface<T,N>::fftwInterface(void) : fftwPlan(NULL) {
}

// destructor
template<typename T, int N>
fftwInterface<T,N>::~fftwInterface(
){
    if ( fftwPlan != NULL ) {
	T tmp = 0;
	destroyPlan(tmp);
    }
}

// create a complex to complex plan;
template <typename T, int N>
int fftwInterface<T,N>::plan( 
    const Array<std::complex<T>, N> &in, 
    Array<std::complex<T>, N> &out, 
    int sign,
    unsigned flags
){
    if ( fftwPlan != NULL ) {
	T tmp = 0;
	destroyPlan(tmp);
    }
    for ( int i = 0; i < N; i++ ) {
        if ( in.length(i) != out.length(i) ) {
	    return -1;
	}
    }
    int* n = new int[N];
    for ( int i = 0; i < N; i++ ) {
	n[i] = in.length(i);
    }
    fftwPlan = planC2C(
		N, 
		n,
		(complex<T>*)in.data(),
		out.data(),
		sign, 
		flags
    );
    return 0;
}

// create a real to complex plan;
template<typename T, int N>
int fftwInterface<T,N>::plan( 
    const Array<T, N> &in, 
    Array<std::complex<T>, N> &out, 
    unsigned flags
) {
    if ( fftwPlan != NULL ) {
	T tmp = 0;
	destroyPlan(tmp);
    }
    int* n = new int[N];
    for ( int i = 0; i < N; i++ ) {
	n[i] = in.length(i);
    }
    fftwPlan = planR2C(
		N, 
		n,
		(T*)in.data(),
		out.data(),
		flags
    );
    return 0;
}

// create a complex to real plan;
template<typename T, int N>
int fftwInterface<T,N>::plan( 
    const Array<std::complex<T>, N> &in, 
    Array<T, N> &out, 
    unsigned flags
) {
    if ( fftwPlan != NULL ) {
	T tmp = 0;
	destroyPlan(tmp);
    }
    int* n = new int[N];
    for ( int i = 0; i < N; i++ ) {
	n[i] = out.length(i);
    }
    fftwPlan = planC2R(
		N, 
		n,
		(complex<T>*)in.data(),
		out.data(),
		flags
    );
    return 0;
}

// create a real to real plan;
template<typename T, int N>
int fftwInterface<T,N>::plan( 
    const Array<T, N> &in, 
    Array<T, N> &out, 
    kindType* kind,
    unsigned flags
) {
    // not implemented yet
	return -1;
}

// execute the plan
template<typename T, int N>
int fftwInterface<T,N>::execute(
){
    if ( fftwPlan != NULL ) {
	T tmp = 0;
        executePlan(tmp);
	return 0;
    } 
    return -1;
}

}
