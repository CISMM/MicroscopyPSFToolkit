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

#ifndef _CONVOLUTION_H
#define _CONVOLUTION_H

#include <blitz/array.h>

using namespace blitz;
namespace cosm {

template<typename T, int N>
Array<T,N> convolution(
    const Array<T,N>& A,
    const Array<T,N>& B
);

template<typename T, int N>
Array<T,N> convolutionCircular(
    const Array<T,N>& A,
    const Array<T,N>& B
);

template<typename T, int N>
Array<T,N> convolutionWithFFT(
	const Array<T,N>& A,
	const Array<T,N>& B,
	unsigned short size = 0,
	//Version 0.8.11::External-3
	//Support an optional check for user 
	//to choose whether Psf is centered.
	bool centeredPsf = true
);

}

#include "convolution.c"

#endif // _CONVOLUTION_H
