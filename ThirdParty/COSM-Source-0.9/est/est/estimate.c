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

#include "estimate.h"

namespace cosm {

// Function to compute an error estimate between two arrays
template<typename T, int N>
T ErrorEstimate<T,N>::error( 
    ErrorType type,
    Array<T,N>& A,
    Array<T,N>& B
) {
    type_ = type;
    return error(A,B);
}

// Function to compute an error estimate between two arrays
template<typename T, int N>
T ErrorEstimate<T,N>::error( 
    Array<T,N>& A,
    Array<T,N>& B
) {
    switch( type_ ) {
      case MAXIMUM_ERROR: return (max)( abs(A - B) );
      case MEAN_ERROR: return (sum)( abs(A - B) )/A.size();
      case MEAN_SQUARE_ERROR: return (sum)( pow2(A-B) )/A.size();
      case LOG_LIKELYHOOD_ERROR:  return logLikelyhood(A,B);
      case I_DIVERGENCE_ERROR: return iDivergence(A,B);
      default:
	return -1;
    }
}

template<typename T, int N>
T ErrorEstimate<T,N>::logLikelyhood( 
    Array<T,N>& A,
    Array<T,N>& B
){
    T epsilon = T(1e-6);
    Array<T,N> C(A.shape());
    C = where ( B < epsilon, epsilon, B );
    return (sum)( B-A * log(C) );
}

template<typename T, int N>
T ErrorEstimate<T,N>::iDivergence( 
    Array<T,N>& A,
    Array<T,N>& B
){
    T epsilon = T(1e-6);
    Array<T,N> C(A.shape());
    Array<T,N> D(A.shape());
    C = where ( A < epsilon, epsilon, A );
    D  = where ( B < epsilon, A, B*(log(B)-log(C))+A-B );
    return (sum)(D);
}
}
