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

#include "blitz/RectDomainIter.h"
#include "blitz/fftwInterface.h"
#include "blitz/arrayManip.h"
#include <blitz/tinyvec-et.h>

namespace cosm {

template<typename T, int N>
Array<T,N> convolution(
    const Array<T,N>& A,
    const Array<T,N>& B
){
    RectDomain<N> domain(A.lbound()+B.lbound(), A.ubound()+B.ubound());
    Array<T,N> C(domain.lbound(), domain.ubound());
    TinyVector<int,N> n(0);
    TinyVector<int,N> m(0);
    RectDomainIter<N> iter(domain);
    for ( n = iter.begin(); !equal(n, iter.end()); n = iter.next() )
    {
        T sum = 0;
	TinyVector<int,N> lb = n - A.ubound();
	TinyVector<int,N> ub = n - A.lbound();
	TinyVector<int,N> low = max(lb, B.lbound());
	TinyVector<int,N> high = min(ub, B.ubound());
    	RectDomain<N> sumDomain(low, high);
        RectDomainIter<N> sumIter(sumDomain);
	for ( m = sumIter.begin(); !equal(m, sumIter.end()); m = sumIter.next() )
	{
	    TinyVector<int,N> k = n - m;
	    sum += A(k)*B(m);
	}
	C(n) = sum;	
    }
    return C;
}

template<typename T, int N>
Array<T,N> convolutionCircular(
    const Array<T,N>& A,
    const Array<T,N>& B
){
}

template<typename T, int N>
Array<T,N> convolutionWithFFT(
    const Array<T,N>& A,
    const Array<T,N>& B,
	unsigned short size,
	//Version 0.8.11::External-3
	//Support an optional check for user 
	//to choose whether Psf is centered.
	bool centeredPsf
){
    TinyVector<int,N> a = A.length();
    TinyVector<int,N> b = B.length();
    Array<T,N> padA(a+b);
    Array<T,N> padB(a+b);

	padAround(A,padA);
	//Version 0.8.11::External-3
    if (centeredPsf){
		padAround(B,padB);
    }else{
		Array<T,N> tempB = B;
		padCenter(tempB,padB);
    }

	//Version 0.8.11::External-3
    Array<std::complex<T>,N> Afft = forwardFFT(padA);
    Array<std::complex<T>,N> Afft1 = forwardFFT(padA);
    Array<std::complex<T>,N> Afft2 = forwardFFT(padA);
    Array<std::complex<T>,N> Bfft = forwardFFT(padB);
    Afft1 = Afft * Bfft;
    Afft2 = Afft * Bfft;

    Array<T,N> C = inverseFFT(Afft1);

    TinyVector<int, N> c = (a+b+1)/2;

	Array<T,N> D1 = circularShift(C,c);

    Array<T,N> D2 = inverseFFT(Afft2);


	TinyVector<int, N> shift = 0;
	if ( size == 1 ) 
	{
	    shift = c-a/2;
		C.resize(A.shape());
    }
	else if ( size == 2 ) 
	{
	    shift = c-b/2;
		C.resize(B.shape());
	}

	//Version 0.8.11::External-3
	if (centeredPsf){
		transform(D1, C, shift);
	}else{
		transform(D2, C, shift);
	}
	return C;
}

}
