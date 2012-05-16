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

#ifndef _ESTIMATE_H
#define _ESTIMATE_H

#include <blitz/array.h>

using namespace blitz;
namespace cosm {

enum ErrorType {
   UNKNOWN_ERROR 	= 0x00,
   MAXIMUM_ERROR 	= 0x01,   // max( abs(est - old) )
   MEAN_ERROR 		= 0x02,   // sum( abs(est - old) )/size 
   MEAN_SQUARE_ERROR 	= 0x04,   // sum( abs(pow2(est)-pow2(old)) )
   LOG_LIKELYHOOD_ERROR = 0x08,   // log likelyhood error estimate 
   I_DIVERGENCE_ERROR  	= 0x10    // i-divergence error estimate
};

template<typename T, int N>
class ErrorEstimate {

  public:


    // ErrorEstmate class constructor
    ErrorEstimate( ErrorType type = UNKNOWN_ERROR ) : type_(type), error_(-1) {};

    // ErrorEstimate class destructor
    ~ErrorEstimate() {};

    // Function to compute an error estimate against a given image
    T error( 
	Array<T,N>& A,
	Array<T,N>& B 
    );
    // Function to compute an error estimate against a given image
    T error( 
	ErrorType type,
	Array<T,N>& A,
	Array<T,N>& B 
    );

    // Function to get last computed error value
    T value() { return error_; };

    // Function to get last computed error type
    ErrorEstimate<T,N> type() { return type_; };

  private:

    T logLikelyhood(           
        Array<T,N>& A,
        Array<T,N>& B
    );

    T iDivergence( 
	Array<T,N>& A,
	Array<T,N>& B 
    );

  private:

    ErrorType type_;
    T error_;
    
};
    

template<typename T, int N>
class Estimate {

  public:

  
    // Estmate class constructor
    Estimate(
	Array<T,N>& psf,
	Array<T,N>& img
    ) : psf_(psf), img_(img) { 
	est_.resize(img_.extent()); 
    };

    // Estmate class constructor
    Estimate(
	Array<T,N>& img
    ) : img_(img) { 
	est_.resize(img_.extent()); 
    };

    // Estimate class destructor
    virtual ~Estimate() { };

    // Function to run the em Algorithm for a number of iterations
    virtual int run() = 0;

    // Function to specify the psf
    virtual void setPSF( 
	Array<T,N>& psf 
    ) { 
	psf_.resize( psf.extent() ); 
	psf_ = psf;
    };

    // Function to specify the raw image
    virtual void setImage( 
	Array<T,N>& img 
    ) { 
	img_.resize( img.extent() ); 
	img_ = img; 
    };

    // Function to get the current estimate
    Array<T,N>& results( 
    ) { return est_; };

  protected:
    
    Array<T,N> psf_;			// point spread function array
    Array<T,N> img_;			// raw image array
    Array<T,N> est_;			// estimation result array

};

}

#include "estimate.c"


#endif // _ESTIMATE_H
