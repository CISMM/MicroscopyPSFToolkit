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

#ifndef ARRAY_MANIP_H
#define ARRAY_MANIP_H

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include "blitz/RectDomainIter.h"

using namespace blitz;

namespace cosm {

template<typename T, int N>
Array<T,N> normalShift(
    Array<T,N>& A,
    TinyVector<int,N>& shift
){
    TinyVector<int, N> pos;
    TinyVector<int, N> size;
    Array<T,N> B(A.shape());
    //Array< T,N >::iterator iter = A.begin(), end = A.end();
    ArrayIterator<T,N> iter = A.begin(), end = A.end();
    while ( iter != end ) 
    {
        bool flag = true;
	pos = (iter.position() + shift);
        size = A.length();
        for ( int i = 0; i < A.dimensions(); i++ )
        {
            if ( pos(i) > size(i) )
            {
                flag = false;
                break;
            } 
        }
	if ( flag )
        {
	    B(pos) = (*iter);
        }
	++iter;
    } 
    return B;
};

template<typename T, int N>
Array<T,N> circularShift(
    Array<T,N>& A,
    TinyVector<int,N>& shift
){
    TinyVector<int, N> pos;
    Array<T,N> B(A.shape());
    //Array< T,N >::iterator iter = A.begin(), end = A.end();
    ArrayIterator<T,N> iter = A.begin(), end = A.end();
    while ( iter != end ) 
    {
	pos = (iter.position() + shift) % A.length();
	B(pos) = (*iter);
	++iter;
    } 
    return B;
};

template<typename T, int N>
void circularShiftAndPadCenter(
    Array<T,N>& A,
    Array<T,N>& B
){
    TinyVector<int,N> l = A.lbound(); 
    TinyVector<int,N> u = A.ubound() + l;
    TinyVector<int,N> shift = ((B.length() - A.length()+1)/2);
    RectDomain<N> rect(l,u);
    B = 0;
    B(rect) = A;
    B = circularShift(B, shift);
};

template<typename T, int N>
void padCenter(
    Array<T,N>& A,
    Array<T,N>& B
){
    TinyVector<int,N> l = A.lbound(); 
    TinyVector<int,N> u = A.ubound();
    TinyVector<int,N> aShift = (A.length()+1)/2;
    TinyVector<int,N> shift = B.length() - aShift;
    RectDomain<N> rect(l,u);
    B = 0;
    B(rect) = circularShift(A, aShift);
    B = circularShift(B, shift);
};

template<typename T, int N>
void transform(
    Array<T,N>& A,
    Array<T,N>& B,
    TinyVector<int,N> shift,
    T value = 0
){
    RectDomain<N> rectA(
	max(A.lbound(), A.lbound()+shift), 
        min(A.ubound(), A.lbound()+shift+(B.ubound()-B.lbound()))
    );
    RectDomain<N> rectB(
        max(B.lbound(), B.lbound()-shift), 
        min(B.ubound(), B.lbound()-shift+(A.ubound()-A.lbound()))
    );
    B = value;
    B(rectB) = A(rectA);
};

template<typename T, int N>
void downsample(
    Array<T,N>& A,
    Array<T,N>& B,
    TinyVector<int,N> factor
){
    // make the B array correct size
    B.resize(A.extent()/factor);

    // iterate over all values of B
    TinyVector<int,N> pos;
    ArrayIterator<T,N> iter = B.begin(), end = B.end();
    while ( iter != end ) 
    {
        pos = iter.position();
        // assign B the average of A values
        RectDomain<N> rect(pos*factor, (pos+1)*factor-1);
        B(pos) = mean( A(rect) );
	++iter;
    } 
};

template<typename T, int N>
void upsample(
    Array<T,N>& A,
    Array<T,N>& B,
    TinyVector<int,N> factor
){
    // make the B array correct size
    B.resize(A.extent()*factor);

    // iterate over all values of A 
    TinyVector<int,N> pos;
    TinyVector<int,N> sub;
    ArrayIterator<T,N> iter = A.begin(), end = A.end();
    while ( iter != end ) 
    {
        pos = iter.position();
        // Replicate new values in B with same A value
        RectDomain<N> rect(pos*factor, (pos+1)*factor-1);
	B(rect) = A(pos);
	++iter;
    } 
};

template<typename T, int N>
void padAround(
    const Array<T,N>& A,
    Array<T,N>& B
){
    TinyVector<int,N> l = A.lbound(); 
    TinyVector<int,N> u = A.ubound() + l;
    TinyVector<int,N> shift = B.length()/2 - A.length()/2;
    RectDomain<N> rect(l,u);
    B = 0;
    B(rect) = A;
    B = circularShift(B, shift);
};

template<typename T, int N>
void mirror(
    Array<T,N>& A,
    Array<T,N>& B
) {
    //Array<T,N>::iterator iter = A.begin(), end = A.end();
    ArrayIterator<T,N> iter = A.begin(), end = A.end();
    while ( iter != end ) 
    {
	TinyVector<int,N> pos1 = iter.position();
	TinyVector<int,N> pos2 = pos1;
 	pos2(0) = B.length(0) - pos1(0) - 1;
	B(pos1) = (*iter);
	B(pos2) = (*iter);
	++iter;
    } 
};

template<typename T, int N>
void multiplyStratum(
    RectDomain<N> rect,
    Array<T,N>& A,
    Array<T,1>& a,
    bool flag = true
) {
    TinyVector<int,N> n(0);
    RectDomainIter<N> iter(rect);
    for ( n = iter.begin(); !equal(n, iter.end()); n = iter.next() )
    {
	A(n) *= ( flag ) ? a(n(0)) : 1 - a(n(0));	
    }
};


template<typename T, int N>
void addBox(
    Array<T,N>& A,
    Array<T,N>& B,
    TinyVector<int,N> center,
    TinyVector<int,N> halfSize,
    T value
){
    B.resize(A.extent());
	B = A;
    RectDomain<N> rect( min(center-halfSize,A.lbound()), 
                        max(center+halfSize, A.ubound()) );
    B(rect) = value; 
};

template<typename T, int N>
void addEllipsoid(
    Array<T,N>& A,
    Array<T,N>& B,
    TinyVector<int,N> center,
    TinyVector<int,N> radius,
    T value
){
    B.resize(A.extent());
	B = A;
    TinyVector<int,N> xNorm(0,0,1);
    TinyVector<int,N> yNorm(0,1,0);
    TinyVector<int,N> zNorm(1,0,0);
    ArrayIterator<T,N> iter = B.begin(), end = B.end();
    while ( iter != end )
    {
        TinyVector<int,N> pos = iter.position();
        TinyVector<double,N> diff = pos - center;
        double sum = 0;
        for ( int i = 0; i < A.dimensions(); i++ )
        {
            sum += pow2(diff(i)/radius(i));
        }
        if ( sum <= 1.0 )
        {
            B(pos) = value;
        }
        ++iter;
    }
};

};
#endif // ARRAY_MAINIP_H
